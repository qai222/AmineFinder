import logging
from io import StringIO

from mfinder import MoleculeFinder, AmineIndicators, MetalIndicators, AmineSubRules
from opsin import OpsinConverter

OC = OpsinConverter("./opsin-2.5.0.jar", rm_tmp=True, wdir="./", silent=True)

AmineFinder = MoleculeFinder(opsin_converter=OC, findername="AmineFinder", kw_indicators=AmineIndicators,
                             kw_exclusion=MetalIndicators, subrules=AmineSubRules)

if __name__ == '__main__':
    log_stream = StringIO()
    logging.basicConfig(stream=log_stream, level=logging.INFO)

    test_names = [
        "catena-[1,2-Ethyldiamine bis(\u03BC4-phosphato)-di-zinc]",
        "catena-(1,8-diaminooctane bis(bis(\u03BC4-phosphato)-(\u03BC2-hydroxo)-tri-tin))",
    ]

    for n in test_names:
        mol, smiles = AmineFinder.get_molecule_and_smiles(n, google_suggestion=False)
        log = log_stream.getvalue()
        print(log)
        print("final:", mol, smiles)
        print("-" * 9)
