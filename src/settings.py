import logging
import colorlog
from pathlib import Path

# ---------------------------------------- Output Setting----------------------------------------------------------------------
HERE = Path(__file__).parent.absolute()
PlotStyle = HERE / 'visualization/paper.mplstyle'
FigureDir = HERE.parent / 'report/figure'
if not FigureDir.is_dir():
    FigureDir.mkdir()


def init_logger(dunder_name, testing_mode=False) -> logging.Logger:
    log_format = (
        '%(asctime)s - '
        '%(name)s - '
        '%(funcName)s - '
        '%(levelname)s - '
        '%(message)s'
    )
    bold_seq = '\033[1m'
    colorlog_format = f'{bold_seq} ' '%(log_color)s ' f'{log_format}'
    colorlog.basicConfig(format=colorlog_format)
    logger = logging.getLogger(dunder_name)

    if testing_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    # # Output full log
    # fh = logging.FileHandler('app.log')
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    # # Output warning log
    # fh = logging.FileHandler('app.warning.log')
    # fh.setLevel(logging.WARNING)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    # # Output error log
    # fh = logging.FileHandler('app.error.log')
    # fh.setLevel(logging.ERROR)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    return logger


# ---------------------------------------- Palette Setting----------------------------------------------------------------------
import palettable

VIVID_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
COLOR_PAlETTE = {
    'Treatment_Arm': {"Chemo->ICI": '#003049', "ICI->Chemo": "#fca311"},
    'RCB': {
        0: '#006d77',
        '0': '#006d77',
        'I': '#83c5be',
        'II': '#ffddd2',
        'III': '#e29578',
    },
    'BestResponse': {'0-I': '#006d77', 'II-III': '#e29578'},
    'Timepoint': {
        'Baseline': "#8ecae6",
        "W3D1": "#219ebc",
        "W7D1": "#023047",
        "Surg": "#ffb703",
        "Surg+AC": "black",
        "EOT?": "gray",
    },
    'Mutation type': {
        'Missense': VIVID_10[5],
        'Nonsense': VIVID_10[0],
        'In frame indel': VIVID_10[1],
        'Frameshift indel': VIVID_10[4],
        'Splice site': VIVID_10[9],
    },
    'WES_APOBEC_Enriched': {'yes': 'salmon', 'no': 'skyblue'},
    'WES_absolute_purity':{'color':'Purples'},
    'TMB': {'WES_Nonsynonymous_TMB': 'purple', 'WES_Synonymous_TMB': 'pink'},
    'Mutation clonality': {'WES_Clonal_TMB': 'purple', 'WES_Subclonal_TMB': 'pink'},
    "Mutation_Signature": {
        'SBS5': '#d8e2dc',
        'SBS46':'tan',
        'SBS1': '#b8bedd',
        'SBS39': '#f0a6ca',
        'SBS2':'skyblue',
        'SBS15':'wheat',
        'SBS3':'violet',
        'SBS42':'gray',
        'SBS6':'moccasin',
        'SBS96':'lightgray',
    },
    "WES_facets_wgd_bool": {True: 'black', False: 'white'},
    "bulkRNA_PAM50": {
        'Normal': 'steelblue',
        'Basal': 'orange',
        'LumB': 'yellow',
        'LumA': 'lightgreen',
        'Her2': 'purple',
    },
    "BluePrint": {
        'Normal': 'steelblue',
        'Basal': 'orange',
        'LumB': 'yellow',
        'LumA': 'lightgreen',
        'Her2': 'purple',
    },
    "er_status":{
        "Weakly Positive (1-10% cell staining)":"red",
        "Positive (>10% cell staining)":"green"
    },
    
}

## Update N/A category
for _, v in COLOR_PAlETTE.items():
    v.update({'N/A': '#e5e5e5'})
