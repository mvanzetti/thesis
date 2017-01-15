banner_default_message = '2017, Manuel Vanzetti'

main_directory = '/Users/manuel/development/thesis/'

etl_phase = [{
    'download': 'download',
    'staging': 'staging',
    'overlap': 'overlap',
    'export': 'export'
}]

encyclopedia = [{
    'encode': {
        'name': 'ENCODE',
        'version': '3',
        'url': 'https://www.encodeproject.org/',
    },
    'fantom': {
        'name': 'FANTOM',
        'version': '5',
        'url': 'http://fantom.gsc.riken.jp/5/',
        'data_url': {
            'presets': 'http://enhancer.binf.ku.dk/presets/'
        }
    },
    'dbsuper': {
        'name': 'dbSUPER',
        'version': '1',
        'url': 'http://bioinfo.au.tsinghua.edu.cn/dbsuper/'
    },
    'epigenomics_roadmap': {

    }
}]
