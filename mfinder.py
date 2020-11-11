import json
import logging
import os
import re
import urllib
from urllib.parse import quote

import yaml
from Levenshtein import ratio

from opsin import OpsinConverter

"""
let <name> be a long IUPAC name contains one or more amine names (<AmineNames>)

1. tokenize <name> to find the tokens (<hits>) containing any of the amine-indicating keywords (<AmineIndicators>)
2. extend these tokens towards left and/or right to get a list of possible amine-containing words, it does NOT 
    stop extending even if it hits a space
3. try to convert all amine-containing words to smiles using opsin, exclude failed cases

Such procedures can be applied to other types of molecules, as long as a set of keywords similar to <AmineIndicators>
can be defined.

google correction: 
    in step 3. there is an option to use google suggestion API to correct possible typos present in <name>, which,
    unfortunately, is not uncommon in the Cambridge Structure Database. 
    The problem is google only tells you the most popular alternatives. You should always check the results.
    Improvements could be expected if other chemical-related APIs (pubchem, chemspider, etc.) are used.

namesub:
    one can define a set of substitution rules to be applied on the <name> before step 1., e.g. `csdsub.yml`
    

"""
logger = logging.getLogger(__name__)
this_folder = os.path.dirname(os.path.abspath(__file__))

BRAS = ["(", "[", "{"]
KETS = [")", "]", "}"]

AmineIndicators = ("piperaz", "benzidin", "eniminium", "adeninium", "lidinium", "linium", "naphthyridin",
                   "phenanthrolin", "piperid", "guanid", "amin", "ammon", "ine", "ammin", "az", "pyrimid", "nitrilo",
                   "phyrin", "pyrrol", "isoindole", "anil", "ptilocaulin", "phycene", "annulene", "metformin", "pyrid")

MetalIndicators = (
    "Actinium", "Aluminium", "Americium", "Barium", "Berkelium", "Beryllium", "Bismuth", "Bohrium", "Cadmium",
    "Calcium", "Californium", "Cerium", "Cesium", "Chromium", "Cobalt", "Copper", "Curium", "Darmstadtium", "Dubnium",
    "Dysprosium", "Einsteinium", "Erbium", "Europium", "Fermium", "Francium", "Gadolinium", "Gallium", "Germanium",
    "Gismondine", "Gold", "Hafnium", "Hassium", "Holmium", "Hydronium", "Hydroxonium", "Indium", "Iridium", "Iron",
    "Lanthanum", "Lawrencium", "Lead", "Lithium", "Lutetium", "Magnesium", "Manganese", "Meitnerium", "Mendelevium",
    "Mercury", "Molybdenum", "Neodymium", "Neptunium", "Nickel", "Niobium", "Nobelium", "Osmium", "Oxonium",
    "Palladium", "Platinum", "Plutonium", "Polonium", "Potassium", "Praseodymium", "Promethium", "Protactinium",
    "Radium", "Rhenium", "Rhodium", "Roentgenium", "Rubidium", "Ruthenium", "Rutherfordium", "Samarium", "Scandium",
    "Seaborgium", "Silver", "Sodium", "Strontium", "Tantalum", "Technetium", "Tellurium", "Terbium", "Thallium",
    "Thorium", "Thulium", "Tin", "Titanium", "Tungsten", "Ununbium", "Ununhexium", "Ununpentium", "Ununquadium",
    "Ununtrium", "Uranium", "Vanadium", "Ytterbium", "Yttrium", "Zinc", "Zirconium"
)

with open("{}/csdsub.yml".format(this_folder), 'r') as stream:
    AmineSubRules = yaml.safe_load(stream)


def StringHits(s: str, ks: [str], excluding=None, regexformat=r"[a-zA-Z0-9]*{}[a-zA-Z0-9]*"):
    """
    find tokens in which there is at least one of <ks>
    :param excluding: a list of keywords that should not appear in the tokens
    """
    ms = []
    excluding = [ex.lower() for ex in excluding]
    for k in ks:
        regex = regexformat.format(k)
        ms += re.findall(regex, s, re.IGNORECASE)
    if excluding is None:
        return set(ms)
    else:
        return set([hit for hit in ms if hit.lower() not in excluding])


def ReplaceByDict(m: str, d: dict, case_insensitive=True):
    nm = m
    for old, new in d.items():
        if case_insensitive:
            iregex = re.compile(re.escape(old), re.IGNORECASE)
            nm = iregex.sub(new, nm)
        else:
            nm = nm.replace(old, new)
    return nm


def GoogleSuggestion(name: str):
    """use google suggestion and Levenshtein ratio to find the correction, return None if failed"""
    qname = "\"{}\"".format(name)
    url = "http://suggestqueries.google.com/complete/search?output=firefox&hl=en&q=" + quote(qname)
    req = urllib.request.Request(url)
    try:
        # req = urllib.parse.quote(req)
        resp = urllib.request.urlopen(req)
    except urllib.error.URLError as e:
        logger.exception(e)
        resp = None
    try:
        resp = resp.read().decode(resp.headers.get_content_charset() or 'utf-8')
        resplist = json.loads(resp)[1]
        result = []
        for s in resplist:
            result += re.findall(r"\w+", s)
        # result = [n for n in json.loads(resp)[1] if " " not in n.strip() and n.strip() != name]
        if len(result) == 0:
            return None
        else:
            result = sorted(result, key=lambda x: ratio(qname, x), reverse=True)
            return result[0]
    except Exception as e:
        logger.exception(e)
        return None


class MoleculeFinderError(Exception): pass


class MoleculeFinder:

    def __init__(self, opsin_converter: OpsinConverter, findername: str, kw_indicators: tuple or list,
                 kw_exclusion: tuple or list = None, subrules: dict = None):
        self.opsin_converter = opsin_converter
        self.findername = findername
        self.kw_indicators = kw_indicators
        self.kw_exclusion = kw_exclusion
        if subrules is None:
            self.subrules = dict()
        else:
            self.subrules = subrules

    def __repr__(self):
        return self.findername

    @staticmethod
    def clean(name: str):
        # avoid additional spaces,
        # e.g. DADTUY12 : catena-[dimethylammonium tris(μ-formato)-cobalt(ii) ]
        for bk in BRAS + KETS:
            p = r"\s*\{}\s*".format(bk)
            name = re.sub(p, bk, name)
        return name

    def sub(self, name, case_insensitive=True):
        name = ReplaceByDict(name, self.subrules, case_insensitive)
        logger.warning("name sub: --> \n \t{}".format(name))
        return name

    def get_hits(self, name):
        return sorted(StringHits(name, self.kw_indicators, excluding=self.kw_exclusion))

    @staticmethod
    def get_extended_hit(name, hit: str):
        # there could be another way to fix typo: for each token find alternatives with google suggestion and use them to extend hits
        # problem is this leads to exponentially large results
        ileft = name.find(hit)
        iright = ileft + len(hit)
        left = name[:ileft]
        right = name[iright:]
        # right = self.name[iright:].split()[0]
        ltokens = [t for t in re.split(r'(\W)', left) if len(t) > 0]
        rtokens = [t for t in re.split(r'(\W)', right) if len(t) > 0]
        possible_lstrings = [""]
        substring = ""
        for t in reversed(ltokens):
            substring = t + substring
            possible_lstrings.append(substring)
        possible_rstrings = [""]
        substring = ""
        for t in rtokens:
            substring = substring + t
            possible_rstrings.append(substring)
        extended_hits = [hit]
        for lstr in possible_lstrings:
            for rstr in possible_rstrings:
                extended_hits.append(lstr + hit + rstr)
        return extended_hits

    def get_all_extended_hits(self, name, hits):
        aehs = []
        for h in hits:
            aehs += self.get_extended_hit(name, h)
        return aehs

    def sniff_check(self, name):
        hits = self.get_hits(name)
        if len(hits) == 0:
            emsg = "no possible hit by {}!: {}".format(self.findername, name)
            raise MoleculeFinderError(emsg)
        possible_strings = self.get_all_extended_hits(name, hits)
        if len(possible_strings) == 0:
            emsg = "no possible extended hit by {}!: {}".format(self.findername, name)
            raise MoleculeFinderError(emsg)
        return hits, sorted(set(possible_strings))

    @staticmethod
    def google_correct_name_by_hits(name: str, hits: [str]):
        """for each hit, try to correct it and replace it in the name"""
        corrected_names = []
        for hit in hits:
            corrected_hit = GoogleSuggestion(hit)
            if corrected_hit is not None:
                corrected_names.append(name.replace(hit, corrected_hit))
        return sorted(set(corrected_name for corrected_name in corrected_names if corrected_name != name))

    def parse_one_name(self, name):
        name = self.clean(name)
        name = self.sub(name)
        hits, possible_strings = self.sniff_check(name)
        logger.info("hits:\n\t{}".format(hits))
        name2smiles = self.opsin_converter.convert_names(possible_strings, output="smi")
        logger.info("converted:\n\t{}".format(name2smiles))
        return name2smiles, hits

    def parse_name(self, name, google_correction=False):
        logger.info("*** parse name in {}: {}".format(self.findername, name))
        name2smiles, hits = self.parse_one_name(name)
        if len(name2smiles) == 0:
            if not google_correction:
                emsg = "no IUPAC names by {}!: {}".format(self.findername, name)
                raise MoleculeFinderError(emsg)
            else:
                emsg = "no IUPAC names by {}!: {} \ntrying google correct, this is unreliable!!".format(self.findername,
                                                                                                        name)
                logger.warning(emsg)
                corrected_names = self.google_correct_name_by_hits(name, hits)
                if len(corrected_names) == 0:
                    emsg = "google cannot suggest new name!"
                    raise MoleculeFinderError(emsg)
                corrected_name2smiles = dict()
                for corrected_name in corrected_names:
                    logger.info("parsing corrected name: {}".format(corrected_name))
                    name2smiles_corrected, hits_from_corrected = self.parse_one_name(corrected_name)
                    for k, v in name2smiles_corrected.items():
                        corrected_name2smiles[k] = v
                if len(corrected_name2smiles) == 0:
                    emsg = "even google suggestion cannot save this: {}".format(name)
                    raise MoleculeFinderError(emsg)
                return corrected_name2smiles
        return name2smiles

    def get_molecule_and_smiles(self, name, google_suggestion=False):
        """get the molecule with the longest name"""
        name2smiles = self.parse_name(name, google_suggestion)
        mol_name = ""
        for n, smi in name2smiles.items():
            if "." in smi:
                continue
            if len(n) > len(mol_name):
                mol_name = n
        return mol_name, name2smiles[mol_name]
