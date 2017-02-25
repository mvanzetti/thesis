overlap_tests = [{
    'fisher_jaccard': {
        'columns': ['encyclopedia', 'biosample_name', 'ovlp_encyclopedia', 'encyclopedia_size',
                    'ovlp_encyclopedia_size',
                    'min_ovlp', 'ovlp_count', 'fisher_left_p', 'fisher_right_p', 'fisher_two_p', 'fisher_oddsratio',
                    'jaccard']

    },
    'fisher_jaccard_z': {
        'columns': ['encyclopedia', 'biosample_name', 'ovlp_encyclopedia', 'encyclopedia_size',
                    'ovlp_encyclopedia_size',
                    'min_ovlp', 'ovlp_count', 'z_random', 'z_shuffled', 'fisher_right_p', 'jaccard']
    },
    'reldist': {
        'columns': ['encyclopedia', 'biosample_name', 'ovlp_encyclopedia', 'encyclopedia_size',
                    'ovlp_encyclopedia_size',
                    'reldist', 'ovlp_count', 'ovlp_fraction']
    }

}]
