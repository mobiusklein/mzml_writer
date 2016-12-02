# http://www.peptideatlas.org/tmp/mzML1.1.0.html

import warnings
import operator
import re

from datetime import datetime
from numbers import Number as NumberBase
from itertools import chain
from functools import partial, update_wrapper
from collections import Mapping

from . import controlled_vocabulary
from .utils import ensure_iterable

from lxml import etree


def make_counter(start=1):
    '''
    Create a functor whose only internal piece of data is a mutable container
    with a reference to an integer, `start`. When the functor is called, it returns
    current `int` value of `start` and increments the mutable value by one.

    Parameters
    ----------
    start: int, optional
        The number to start counting from. Defaults to `1`.

    Returns
    -------
    int:
        The next number in the count progression.
    '''
    start = [start]

    def count_up():
        ret_val = start[0]
        start[0] += 1
        return ret_val
    return count_up


def camelize(name):
    parts = name.split("_")
    if len(parts) > 1:
        return ''.join(parts[0] + [part.title() if part != "ref" else "_ref" for part in parts[1:]])
    else:
        return name


def id_maker(type_name, id_number):
    return "%s_%s" % (type_name.upper(), str(id_number))


def sanitize_id(string):
    string = re.sub(r"\s", '_', string)
    string = re.sub(r"\\|/", '', string)
    return string


NO_TRACK = object()


class CountedType(type):
    _cache = {}

    def __new__(cls, name, parents, attrs):
        new_type = type.__new__(cls, name, parents, attrs)
        tag_name = attrs.get("tag_name")
        new_type.counter = staticmethod(make_counter())
        if attrs.get("_track") is NO_TRACK:
            return new_type
        cls._cache[name] = new_type
        if tag_name is not None:
            cls._cache[tag_name] = new_type
        return new_type


class TagBase(object):
    __metaclass__ = CountedType

    type_attrs = {}

    def __init__(self, tag_name=None, text="", **attrs):
        self.tag_name = tag_name or self.tag_name
        _id = attrs.pop('id', None)
        self.attrs = {}
        self.attrs.update(self.type_attrs)
        self.text = text
        self.attrs.update(attrs)
        # When passing through a XMLWriterMixin.element() call, tags may be reconstructed
        # and any set ids will be passed through the attrs dictionary, but the `with_id`
        # flag won't be propagated. `_force_id` preserves this.
        self._force_id = True
        if _id is None:
            self._id_number = self.counter()
            self._id_string = None
            self._force_id = False
        elif isinstance(_id, int):
            self._id_number = _id
            self._id_string = None
        elif isinstance(_id, basestring):
            self._id_number = None
            self._id_string = _id

    def __getattr__(self, key):
        try:
            return self.attrs[key]
        except KeyError:
            try:
                return self.attrs[camelize(key)]
            except KeyError:
                raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, key))

    # Support Mapping Interface
    def __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key)

    def keys(self):
        yield "id"
        for key in self.attrs:
            yield key

    def get(self, name, default):
        return getattr(self, name, default)

    @property
    def id(self):
        if self._id_string is None:
            self._id_string = id_maker(self.tag_name, self._id_number)
        return self._id_string

    @property
    def with_id(self):
        return False

    def element(self, xml_file=None, with_id=False):
        with_id = self.with_id or with_id or self._force_id
        # if self.tag_name == "software" and not with_id:
        #     raise Exception()
        attrs = {k: str(v) for k, v in self.attrs.items() if v is not None}
        if with_id:
            attrs['id'] = self.id
        if xml_file is None:
            return etree.Element(self.tag_name, **attrs)
        else:
            return xml_file.element(self.tag_name, **attrs)

    def write(self, xml_file, with_id=False):
        el = self.element(with_id=with_id)
        xml_file.write(el)

    __call__ = element

    def __repr__(self):
        return "<%s id=\"%s\" %s>" % (self.tag_name, self.id, " ".join("%s=\"%s\"" % (
            k, str(v)) for k, v in self.attrs.items()))

    def __eq__(self, other):
        try:
            return self.attrs == other.attrs
        except AttributeError:
            return False

    def __ne__(self, other):
        try:
            return self.attrs != other.attrs
        except AttributeError:
            return True

    def __hash__(self):
        return hash((self.tag_name, frozenset(self.attrs.items())))


class MzML(TagBase):
    type_attrs = {
        "xmlns": "http://psi.hupo.org/ms/mzml",
        "version": "1.1.0",
        "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
        "xsi:schemaLocation": "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"
    }

    def __init__(self, **attrs):
        attrs.setdefault('creationDate', datetime.utcnow())
        super(MzML, self).__init__("mzML", **attrs)


class CVParam(TagBase):
    tag_name = "cvParam"

    @classmethod
    def param(cls, name, value=None, **attrs):
        if isinstance(name, cls):
            return name
        else:
            if value is None:
                return cls(name=name, **attrs)
            else:
                return cls(name=name, value=value, **attrs)

    def __init__(self, accession=None, name=None, ref=None, value=None, **attrs):
        if ref is not None:
            attrs["cvRef"] = ref
        if accession is not None:
            attrs["accession"] = accession
        if name is not None:
            attrs["name"] = name
        if value is not None:
            attrs['value'] = value
        else:
            attrs['value'] = ''
        super(CVParam, self).__init__(self.tag_name, **attrs)
        self.patch_accession(accession, ref)

    @property
    def value(self):
        return self.attrs.get("value")

    @value.setter
    def value(self, value):
        self.attrs['value'] = value

    @property
    def ref(self):
        return self.attrs['cvRef']

    @property
    def name(self):
        return self.attrs['name']

    @property
    def accession(self):
        return self.attrs['accession']

    def __call__(self, *args, **kwargs):
        self.write(*args, **kwargs)

    def __repr__(self):
        return "<%s %s>" % (self.tag_name, " ".join("%s=\"%s\"" % (
            k, str(v)) for k, v in self.attrs.items()))

    def patch_accession(self, accession, ref):
        if accession is not None:
            if isinstance(accession, int):
                accession = "%s:%d" % (ref, accession)
                self.attrs['accession'] = accession
            else:
                self.attrs['accession'] = accession


class UserParam(CVParam):
    tag_name = "userParam"


class CV(TagBase):
    tag_name = 'cv'

    def __init__(self, id, uri, **kwargs):
        super(CV, self).__init__(id=id, uri=uri, **kwargs)
        self._vocabulary = None

    def load(self, handle=None):
        if handle is None:
            fp = controlled_vocabulary.obo_cache.resolve(self.uri)
            cv = controlled_vocabulary.ControlledVocabulary.from_obo(fp)
        else:
            cv = controlled_vocabulary.ControlledVocabulary.from_obo(handle)
        try:
            cv.id = self.id
        except:
            pass
        return cv

    def __getitem__(self, key):
        if self._vocabulary is None:
            self._vocabulary = self.load()
        return self._vocabulary[key]


def identity(x):
    return x


def _make_tag_type(name, **attrs):
    return type(name, (TagBase,), {"tag_name": name, "type_attrs": attrs})


def _element(_tag_name, *args, **kwargs):
    try:
        eltype = CountedType._cache[_tag_name]
    except KeyError:
        eltype = _make_tag_type(_tag_name)
    return eltype(*args, **kwargs)


def element(xml_file, _tag_name, *args, **kwargs):
    with_id = kwargs.pop("with_id", False)
    if isinstance(_tag_name, basestring):
        el = _element(_tag_name, *args, **kwargs)
    else:
        el = _tag_name
    return el.element(xml_file=xml_file, with_id=with_id)


default_cv_list = [
    _element(
        "cv", id="PSI-MS",
        uri=("http://psidev.cvs.sourceforge.net/*checkout*/"
             "psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"),
        version="3.79.0", fullName="PSI-MS"),
    _element(
        "cv", id="UO",
        uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo",
        fullName="UNIT-ONTOLOGY"),
]


common_units = {
    "parts per million": "UO:0000169",
    "dalton": "UO:0000221"
}


class ChildTrackingMeta(type):
    def __new__(cls, name, parents, attrs):
        if not hasattr(cls, "_cache"):
            cls._cache = dict()
        new_type = type.__new__(cls, name, parents, attrs)
        cls._cache[name] = new_type
        return new_type


class SpecializedContextCache(dict):
    def __init__(self, type_name):
        self.type_name = type_name

    def __getitem__(self, key):
        try:
            item = dict.__getitem__(self, key)
            return item
        except KeyError:
            if key is None:
                return None
            warnings.warn("No reference was found for %r in %s" % (key, self.type_name), stacklevel=3)
            new_value = id_maker(self.type_name, key)
            self[key] = new_value
            return new_value

    def __repr__(self):
        return '%s\n%s' % (self.type_name, dict.__repr__(self))


class VocabularyResolver(object):
    def __init__(self, vocabularies=None):
        if vocabularies is None:
            vocabularies = default_cv_list
        self.vocabularies = vocabularies

    def get_vocabulary(self, id):
        for vocab in self.vocabularies:
            if vocab.id == id:
                return vocab
        raise KeyError(id)

    def param(self, name, value=None, cv_ref=None, **kwargs):
        accession = kwargs.get("accession")
        if isinstance(name, CVParam):
            return name

        if isinstance(name, Mapping):
            mapping = name
            value = value or mapping.get('value')
            accession = accession or mapping.get("accession")
            cv_ref = cv_ref or mapping.get("cv_ref") or mapping.get("cvRef")
            name = mapping.get('name')

            kwargs.update({k: v for k, v in mapping.items() if k not in ("name", "value", "accession")})

        if cv_ref is None:
            for cv in self.vocabularies:
                try:
                    term = cv[name]
                    name = term["name"]
                    accession = term["id"]
                    cv_ref = cv.id
                except:
                    pass
        if cv_ref is None:
            return UserParam(name=name, value=value, **kwargs)
        else:
            kwargs.setdefault("ref", cv_ref)
            kwargs.setdefault("accession", accession)
            return CVParam(name=name, value=value, **kwargs)

    def term(self, name, include_source=False):
        for cv in self.vocabularies:
            try:
                term = cv[name]
                if include_source:
                    return term, cv
                else:
                    return term
            except:
                pass
        else:
            raise KeyError(name)


class DocumentContext(dict, VocabularyResolver):
    def __init__(self, vocabularies=None):
        dict.__init__(self)
        VocabularyResolver.__init__(self, vocabularies)

    def __missing__(self, key):
        self[key] = SpecializedContextCache(key)
        return self[key]


NullMap = DocumentContext()


class ReprBorrowingPartial(partial):
    """
    Create a partial instance that uses the wrapped callable's
    `__repr__` method instead of a generic partial
    """
    def __init__(self, func, *args, **kwargs):
        self._func = func
        super(ReprBorrowingPartial, self).__init__(func, *args, **kwargs)
        update_wrapper(self, func)

    def __repr__(self):
        return repr(self.func)

    def __getattr__(self, name):
        return getattr(self._func, name)


class ComponentDispatcher(object):
    """
    A container for a :class:`DocumentContext` which provides
    an automatically parameterized version of all :class:`ComponentBase`
    types which use this instance's context.

    Attributes
    ----------
    context : :class:`DocumentContext`
        The mapping responsible for managing the global
        state of all created components.
    """
    def __init__(self, context=None, vocabularies=None):
        if context is None:
            context = DocumentContext(vocabularies=vocabularies)
        else:
            if vocabularies is not None:
                context.vocabularies.extend(vocabularies)
        self.context = context

    def __getattr__(self, name):
        """
        Provide access to an automatically parameterized
        version of all :class:`ComponentBase` types which
        use this instance's context.

        Parameters
        ----------
        name : str
            Component Name

        Returns
        -------
        ReprBorrowingPartial
            A partially parameterized instance constructor for
            the :class:`ComponentBase` type requested.
        """
        component = ChildTrackingMeta._cache[name]
        return ReprBorrowingPartial(component, context=self.context)

    def register(self, entity_type, id):
        """
        Pre-declare an entity in the document context. Ensures that
        a reference look up will be satisfied.

        Parameters
        ----------
        entity_type : str
            An entity type, either a tag name or a component name
        id : int
            The unique id number for the thing registered

        Returns
        -------
        str
            The constructed reference id
        """
        value = id_maker(entity_type, id)
        self.context[entity_type][id] = value
        return value

    @property
    def vocabularies(self):
        return self.context.vocabularies

    def param(self, *args, **kwargs):
        return self.context.param(*args, **kwargs)

    def term(self, *args, **kwargs):
        return self.context.term(*args, **kwargs)

    def get_vocabulary(self, *args, **kwargs):
        return self.context.get_vocabulary(*args, **kwargs)

# ------------------------------------------
# Base Component Definitions


class ComponentBase(object):
    __metaclass__ = ChildTrackingMeta

    def __init__(self, *args, **kwargs):
        pass

    def __getattr__(self, key):
        try:
            return self.element.attrs[key]
        except KeyError:
            raise AttributeError(key)

    def write(self, xml_file):
        raise NotImplementedError()

    def __call__(self, xml_file):
        self.write(xml_file)

    def __repr__(self):
        return "%s\n%s" % (
            self.element, "\n".join([
                "  %s: %r" % (k, v) for k, v in self.__dict__.items()
                if k not in ("context", "element") and not k.startswith("_")])
        )


class ParameterContainer(ComponentBase):
    def __init__(self, tag_name, params=None, context=NullMap):
        if params is None:
            params = []
        self.element = _element(tag_name)
        self.context = context
        self.params = params

    def write(self, xml_file):
        with self.element(xml_file, with_id=False):
            for param in self.params:
                self.context.param(param)(xml_file)


class GenericCollection(ComponentBase):
    def __init__(self, tag_name, members, context=NullMap):
        self.members = members
        self.tag_name = tag_name
        self.element = _element(tag_name, xmlns="http://psidev.info/psi/pi/mzML/1.1", count=len(self.members))

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            for member in self.members:
                member.write(xml_file)


class IDMemberGenericCollection(GenericCollection):
    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            for member in self.members:
                member.write(xml_file, with_id=True)


class IDGenericCollection(GenericCollection):
    def __init__(self, tag_name, members, id, context=NullMap):
        self.members = members
        self.tag_name = tag_name
        self.element = _element(tag_name, xmlns="http://psidev.info/psi/pi/mzML/1.1", id=id, count=len(self.members))
        context[tag_name][id] = self.element.id

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for member in self.members:
                member.write(xml_file)


# --------------------------------------------------
# File Metadata


class FileContent(ComponentBase):
    def __init__(self, spectrum_types, context=NullMap):
        self.spectrum_types = spectrum_types
        self.element = _element("fileContent")
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file):
            for spectrum_type in self.spectrum_types:
                self.context.param(spectrum_type)(xml_file)


class SourceFileList(GenericCollection):
    def __init__(self, members, context):
        super(SourceFileList, self).__init__("sourceFileList", members, context)


class SourceFile(ComponentBase):
    def __init__(self, location, name, id=None, params=None, context=NullMap):
        if params is None:
            params = []
        self.location = location
        self.name = name
        self.element = _element("SourceFile", location=location, id=id, name=name)
        self.params = params
        self.context = context
        context["SourceFile"][id] = self.element.id

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)


class FileDescription(ComponentBase):
    def __init__(self, content, source_files, contacts=None, context=NullMap):
        if not isinstance(content, FileContent):
            content = FileContent(content, context=context)
        if not isinstance(source_files, SourceFileList):
            if len(source_files) > 0 and not isinstance(source_files[0], SourceFile):
                source_files = [SourceFile(context=context, **f) for f in source_files]
            source_files = SourceFileList(source_files, context=context)
        self.content = content
        self.source_files = source_files
        self.contacts = contacts
        self.context = context
        self.element = _element("fileDescription")

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            self.content.write(xml_file)
            self.source_files.write(xml_file)
            # TODO: handle contact


# --------------------------------------------------
# ParamGroups

class ReferenceableParamGroupList(IDGenericCollection):
    def __init__(self, members, context):
        super(ReferenceableParamGroupList, self).__init__(
            "referenceableParamGroupList", members, context)


class ReferenceableParamGroup(ComponentBase):
    def __init__(self, params=None, id=None, context=NullMap):
        if params is None:
            params = []
        self.params = params
        self.element = _element("referenceableParamGroup", id=id)
        self.id = self.element.id
        self.context = context
        context["ReferenceableParamGroup"][id] = self.element.id

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)


# --------------------------------------------------
# Sample Metadata


class SampleList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(SampleList, self).__init__("sampleList", members, context)


class Sample(ComponentBase):
    def __init__(self, name, params=None, id=None, context=NullMap):
        if params is None:
            params = []
        self.name = name
        self.params = params
        self.element = _element("sample", name=name, id=id)
        self.id = self.element.id
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)

# --------------------------------------------------
# Software Processing Metadata


class SoftwareList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(SoftwareList, self).__init__("softwareList", members, context)


class Software(ComponentBase):
    def __init__(self, id=None, version="0.0", params=None, context=NullMap):
        if params is None:
            params = []
        self.version = version
        self.params = params
        self.element = _element("software", id=id, version=version)
        self.id = self.element.id
        self.context = context
        context['Software'][id] = self.element.id

    @property
    def with_id(self):
        return True

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)


# --------------------------------------------------
# Scan Settings Metadata


class ScanSettingsList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(ScanSettingsList, self).__init__("scanSettingsList", members, context)


class ScanSettings(ComponentBase):
    def __init__(self, id=None, source_file_references=None, target_list=None, params=None, context=NullMap):
        if source_file_references is None:
            source_file_references = []
        if target_list is None:
            target_list = []
        if params is None:
            params = []
        self.params = params
        self.source_file_references = source_file_references
        self.target_list = target_list
        self.element = _element("scanSettings", id=id)
        self.id = self.element.id
        self.context = context
        context['ScanSettings'][id] = self.id

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            source_refs = GenericCollection(
                "sourceFileRefList",
                [_element("sourceFileRef", ref=i) for i in self.source_file_references])
            if len(source_refs):
                source_refs.write(xml_file)
            # TODO handle targetList and targets
            for param in self.params:
                self.context.param(param)(xml_file)


# --------------------------------------------------
# Instrument Configuration Metadata

class InstrumentConfigurationList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(InstrumentConfigurationList, self).__init__("instrumentConfigurationList", members, context)


class InstrumentConfiguration(ComponentBase):
    def __init__(self, scan_settings_reference, id, component_list=None, params=None,
                 software_reference=None, context=NullMap):
        self.scan_settings_reference = scan_settings_reference
        self.params = params
        self.software_reference = software_reference
        self._software_reference = context['Software'][software_reference]
        self.component_list = component_list
        self.element = _element(
            "instrumentConfiguration", id=id,
            scanSettingsRef=context['ScanSettings'][scan_settings_reference])
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)
            if self.component_list is not None:
                self.component_list.write(xml_file)
            if self.software_reference is not None:
                _element("softwareRef", ref=self._software_reference).write(xml_file)


class ComponentList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(ComponentList, self).__init__("componentList", members, context)

    @classmethod
    def build(cls, members, context=NullMap, type_key='type'):
        components = []
        for component in members:
            if component[type_key] == 'source':
                components.append(
                    Source(component['order'], component['params'], context))
            elif component[type_key] == 'analyzer':
                components.append(
                    Analyzer(component['order'], component['params'], context))
            elif component[type_key] == 'detector':
                components.append(
                    Detector(component['order'], component['params'], context))
            else:
                raise KeyError("Unknown component %s" % component[type_key])
        return cls(components, context=context)


class Source(ComponentBase):
    def __init__(self, order, params=None, context=NullMap):
        self.order = order
        self.params = params
        self.element = _element("source", order=order)
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file):
            for param in self.params:
                self.context.param(param)(xml_file)


class Analyzer(ComponentBase):
    def __init__(self, order, params=None, context=NullMap):
        self.order = order
        self.params = params
        self.element = _element("analyzer", order=order)
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file):
            for param in self.params:
                self.context.param(param)(xml_file)


class Detector(ComponentBase):
    def __init__(self, order, params=None, context=NullMap):
        self.order = order
        self.params = params
        self.element = _element("detector", order=order)
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file):
            for param in self.params:
                self.context.param(param)(xml_file)


# --------------------------------------------------
# Data Processing Metadata


class DataProcessingList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(DataProcessingList, self).__init__("dataProcessingList", members, context)


class DataProcessing(ComponentBase):
    def __init__(self, processing_methods=None, id=None, context=NullMap):
        if processing_methods is None:
            processing_methods = []
        if processing_methods and not isinstance(processing_methods[0], ProcessingMethod):
            processing_methods = [ProcessingMethod(context=context, **m) for m in processing_methods]
        self.processing_methods = processing_methods
        self.element = _element("dataProcessing", id=id)
        self.context = context
        context['DataProcessing'][id] = self.element.id

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for method in self.processing_methods:
                method.write(xml_file)


class ProcessingMethod(ComponentBase):
    def __init__(self, order, software_reference, params=None, context=NullMap):
        if params is None:
            params = []
        self.order = order
        self.software_reference = software_reference
        self._software_reference = context['Software'][software_reference]
        self.element = _element("processingMethod", order=order, softwareRef=self._software_reference)
        self.params = params
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)

# --------------------------------------------------
# Spectral and Chromatographic Data Storage


class SpectrumList(ComponentBase):
    def __init__(self, members, default_data_processing_reference, context=NullMap):
        self.members = members
        self.default_data_processing_reference = default_data_processing_reference
        self._default_data_processing_reference = context["DataProcessing"][default_data_processing_reference]
        self.element = _element(
            "spectrumList", count=len(self.members), defaultDataProcessingRef=self._default_data_processing_reference)
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            for member in self.members:
                member.write(xml_file)


class Spectrum(ComponentBase):
    def __init__(self, index, binary_data_list=None, scan_list=None, precursor_list=None, product_list=None,
                 default_array_length=None, source_file_reference=None, data_processing_reference=None,
                 id=None,
                 params=None, context=NullMap):
        if params is None:
            params = []
        self.index = index
        self.scan_list = scan_list
        self.precursor_list = precursor_list
        self.product_list = product_list
        self.binary_data_list = binary_data_list
        self.default_array_length = default_array_length
        self.source_file_reference = source_file_reference
        self._source_file_reference = context["SourceFile"][source_file_reference]
        self.data_processing_reference = data_processing_reference
        self._data_processing_reference = context["DataProcessing"][data_processing_reference]
        self.element = _element(
            "spectrum", id=id, index=index, sourceFileRef=self._source_file_reference,
            defaultArrayLength=self.default_array_length, dataProcessingRef=self._data_processing_reference)
        self.context = context
        self.context["Spectrum"][id] = self.element.id
        self.params = params

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            for param in self.params:
                self.context.param(param)(xml_file)
            if self.scan_list is not None:
                self.scan_list.write(xml_file)
            if self.precursor_list is not None:
                self.precursor_list.write(xml_file)

            self.binary_data_list.write(xml_file)


class Run(ComponentBase):
    def __init__(self, default_instrument_configuration_reference, spectrum_list=None, chromatogram_list=None, id=None,
                 default_source_file_reference=None, sample_reference=None, start_time_stamp=None, params=None,
                 context=NullMap):
        if params is None:
            params = []
        if spectrum_list is None:
            spectrum_list = []
        if chromatogram_list is None:
            chromatogram_list = []
        self.params = params
        self.spectrum_list = spectrum_list
        self.chromatogram_list = chromatogram_list
        self.sample_reference = sample_reference
        self._sample_reference = context["Sample"][sample_reference]
        self.default_instrument_configuration_reference = default_instrument_configuration_reference
        self._default_instrument_configuration_reference = context["InstrumentConfiguration"][
            default_instrument_configuration_reference]
        self.default_source_file_reference = default_source_file_reference
        self._default_source_file_reference = context["SourceFile"][default_source_file_reference]
        self.start_time_stamp = start_time_stamp
        self.element = _element(
            "run", id=id, defaultInstrumentConfigurationRef=self._default_instrument_configuration_reference,
            defaultSourceFileRef=self._default_source_file_reference, sampleRef=self._sample_reference,
            startTimeStamp=start_time_stamp)
        self.context = context


class BinaryDataArrayList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(BinaryDataArrayList, self).__init__("binaryDataArrayList", members, context)


class BinaryDataArray(ComponentBase):
    def __init__(self, binary, encoded_length, data_processing_reference=None, array_length=None,
                 params=None, context=NullMap):
        if params is None:
            params = []
        self.encoded_length = encoded_length
        self.data_processing_reference = data_processing_reference
        self._data_processing_reference = context["DataProcessing"][data_processing_reference]
        self.array_length = array_length
        self.params = params
        self.binary = binary

        self.element = _element(
            "binaryDataArray", arrayLength=array_length, encodedLength=encoded_length,
            dataProcessingRef=self._data_processing_reference)
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            for param in self.params:
                self.context.param(param)(xml_file)
            self.binary.write(xml_file)


class Binary(ComponentBase):
    def __init__(self, encoded_array, context=NullMap):
        self.encoded_array = encoded_array
        self.context = context
        self.element = _element("binary", text=encoded_array)

    def write(self, xml_file):
        with self.element(xml_file, with_id=False):
            xml_file.write(self.encoded_array)


class ScanList(ComponentBase):
    def __init__(self, members, params=None, context=NullMap):
        if params is None:
            params = []
        self.members = members
        self.params = params
        self.element = _element("scanList", count=len(self.members))
        self.context = context

    def write(self, xml_file):
        with self.element(xml_file, with_id=False):
            for param in self.params:
                self.context.param(param)(xml_file)
            for member in self.members:
                member.write(xml_file)


class Scan(ComponentBase):
    def __init__(self, scan_window_list=None, params=None, context=NullMap):
        if scan_window_list is None:
            scan_window_list = ScanWindowList([], context)
        elif not isinstance(scan_window_list, ScanWindowList):
            scan_window_list = ScanWindowList(scan_window_list, context)
        if params is None:
            params = []
        self.params = params
        self.scan_window_list = scan_window_list
        self.element = _element("scan")
        self.context = context

    def write(self, xml_file):
        with self.element(xml_file, with_id=False):
            for param in self.params:
                self.context.param(param)(xml_file)
            self.scan_window_list.write(xml_file)


class ScanWindowList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(ScanWindowList, self).__init__('scanWindowList', members, context)


class ScanWindow(ParameterContainer):
    def __init__(self, *args, **kwargs):
        super(ScanWindow, self).__init__("scanWindow", *args, **kwargs)


class PrecursorList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(PrecursorList, self).__init__('precursorList', members, context)


class Precursor(ComponentBase):
    def __init__(self, selected_ion_list, activation, isolation_window=None, spectrum_reference=None, context=NullMap):
        self.selected_ion_list = selected_ion_list
        self.activation = activation
        self.isolation_window = isolation_window
        self.spectrum_reference = spectrum_reference
        self._spectrum_reference = context["Spectrum"][spectrum_reference]
        self.element = _element("precursor", spectrumRef=self._spectrum_reference)

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=False):
            if self.activation is not None:
                self.activation.write(xml_file)
            if self.isolation_window is not None:
                self.isolation_window.write(xml_file)
            self.selected_ion_list.write(xml_file)


class SelectedIonList(GenericCollection):
    def __init__(self, members, context=NullMap):
        super(SelectedIonList, self).__init__("selectedIonList", members, context)


class SelectedIon(ComponentBase):
    def __init__(self, selected_ion_mz, intensity=None, charge=None, params=None, context=NullMap):
        if params is None:
            params = []
        self.selected_ion_mz = selected_ion_mz
        self.intensity = intensity
        self.charge = charge
        self.params = params
        self.element = _element("selectedIon")
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file):
            if self.selected_ion_mz is not None:
                self.context.param(name="selected ion m/z", value=self.selected_ion_mz)(xml_file)
            if self.intensity is not None:
                self.context.param(name="peak intensity", value=self.intensity)(xml_file)
            if self.charge is not None:
                self.context.param(name="charge state", value=self.charge)(xml_file)


# --------------------------------------------------
# Misc. Providence Management


DEFAULT_CONTACT_ID = "PERSON_DOC_OWNER"
DEFAULT_ORGANIZATION_ID = "ORG_DOC_OWNER"


class CVList(ComponentBase):
    def __init__(self, cv_list=None, context=NullMap):
        if cv_list is None:
            cv_list = default_cv_list
        self.cv_list = cv_list

    def write(self, xml_file):
        with element(xml_file, 'cvList'):
            for member in self.cv_list:
                xml_file.write(member.element(with_id=True))


class Person(ComponentBase):
    def __init__(self, first_name='first_name', last_name='last_name', id=DEFAULT_CONTACT_ID,
                 affiliation=DEFAULT_ORGANIZATION_ID, context=NullMap):
        self.first_name = first_name
        self.last_name = last_name
        self.id = id
        self.affiliation = affiliation
        self.element = _element("Person", firstName=first_name, last_name=last_name, id=id)
        context["Person"][id] = self.element.id
        self.context = context

    def write(self, xml_file):
        with self.element.element(xml_file, with_id=True):
            element(xml_file, 'Affiliation', organization_ref=self.affiliation)


class Organization(ComponentBase):
    def __init__(self, name="name", id=DEFAULT_ORGANIZATION_ID, context=NullMap):
        self.name = name
        self.id = id
        self.element = _element("Organization", name=name, id=id)
        context["Organization"][id] = self.id
        self.context = context

    def write(self, xml_file):
        xml_file.write(self.element.element())


DEFAULT_PERSON = Person()
DEFAULT_ORGANIZATION = Organization()
