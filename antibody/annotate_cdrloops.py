import os
from anarci import anarci,number
from pprint import pprint

light_sequence="DIELTQSPAIMSASPGEKVTMTCSASSSVSYMHWYQQKSGTSPKRWIYDTSKLASGVPGRFSGSGSGNSYSLTISSVEAEDDATYYCQQWSKHPLTFGSGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE"
heavy_sequence="QVQLQQSGPELEKPGASVKISCKASGYSFTGYTMNWVKQSHGKSLEWIGLITPYNGASSYNQKFRGKATLTVDKSSSTAYMDLLSLTSEDSAVYFCARGGYDGRGFDYWGSGTPVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"

nanobody_sequence="QVQLVESGGGLVQAGGSLRLSCAASGFPVAYKTMWWYRQAPGKEREWVAAIESYGIKWTRYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCIVWVGAQYHGQGTQVTVSA"
notantibody_sequence="MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL"
sequences=[("7UED_L",light_sequence),
           ("7UED_H",heavy_sequence),
           ('102L',notantibody_sequence)]
result=anarci(sequences,scheme='imgt',output=False)
result2=anarci(sequences,scheme='chothia',output=False)
#pprint(result)
def annotate_cdrloops(anarci_result):
    """
    annotate the cdr loops from ANARCI output
    Args:
        anarci_result:

    Returns:
        CDRs (List[dict]):
    """
    length=len(anarci_result[0])
    CDRS=[]
    for i in range(length):
        #print(anarci_result[0][i],anarci_result[1][i],anarci_result[2][i])
        if anarci_result[1][i]:
            chain_type=anarci_result[1][i][0]['chain_type']
            query_name=anarci_result[1][i][0]['query_name']
            scheme=anarci_result[1][i][0]['scheme']
            cdrs = {"CDR1": [],
                    "CDR2": [],
                    "CDR3": [], }
            if scheme=='imgt':
                for (pos,_),res in anarci_result[0][i][0][0]:
                    if pos is None or res == '-':
                        continue
                    if 27<=pos<=38:
                        cdrs["CDR1"].append(res)
                    elif 56<=pos<=65:
                        cdrs["CDR2"].append(res)
                    elif 105<=pos<=117:
                        cdrs["CDR3"].append(res)
            elif scheme=='chothia':
                if chain_type == "H":
                    for (pos,_),res in anarci_result[0][i][0][0]:
                        if pos is None or res == '-':
                            continue
                        if 26<=pos<33:
                            cdrs["CDR1"].append(res)
                        elif 52<=pos<57:
                            cdrs["CDR2"].append(res)
                        elif 95<=pos<103:
                            cdrs["CDR3"].append(res)
                elif chain_type == "L" or "K":
                    for (pos,_),res in anarci_result[0][i][0][0]:
                        if pos is None or res == '-':
                            continue
                        if 24<=pos<35:
                            cdrs["CDR1"].append(res)
                        elif 50<=pos<57:
                            cdrs["CDR2"].append(res)
                        elif 89<=pos<98:
                            cdrs["CDR3"].append(res)
                else:
                    raise ValueError(f"Unknown chain type: {chain_type}")

            else:
                raise TypeError("Unknown scheme type, currently only support imgt or chothia")

            for name, residues in cdrs.items():
                cdrs[name]=''.join(residues)

            CDRS.append({'query_name':query_name,'scheme':scheme,'chain_type':chain_type,
                          'CDRS':cdrs})

        else:
            CDRS.append({'query_name':f"_index_{i}",'scheme':None,'chain_type':None,
                         'CDRS':None})

    return CDRS


print(annotate_cdrloops(result2))

