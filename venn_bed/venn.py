#!/usr/bin/env python
from __future__ import print_function

import re
import subprocess
import json
import itertools
import six
from collections import defaultdict
from os.path import dirname, basename, join, splitext, isfile, abspath

from ngs_utils import logger
from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file
from ngs_utils.reporting.reporting import write_static_html_report


def find_pairwise_intersections(work_dirpath, fpaths, regions_file=None):
    if all(fp.endswith('.vcf.gz') for fp in fpaths):
        # For all bgzipped VCFs, we can use effective bcftools isec tool
        intersection_size_by_subset = _intersect_vcfs(work_dirpath, fpaths, regions_file=regions_file)

    else:
        if regions_file:
            fpaths.append(regions_file)
        intersection_file_by_subset = dict()
        _calc_set_intersection(work_dirpath, sorted(fpaths), intersection_file_by_subset)
        intersection_size_by_subset = dict()
        # label_by_subset = dict()
        for file_set, intersection_file in intersection_file_by_subset.items():
            file_set = tuple([splitext(basename(b))[0] for b in file_set])
            if intersection_file.endswith('.bed'):
                intersection_size_by_subset[file_set] = bedsize(intersection_file)
            else:
                intersection_size_by_subset[file_set] = vcfsize(intersection_file)
            # label_by_subset[bed_set] = basename(splitext(intersection_bed)[0])
            print(str(file_set) + ': ' + basename(intersection_file) + ', size: ' + str(intersection_size_by_subset[file_set]))

    return intersection_size_by_subset


def get_upset_file(work_dir, output_dir, fpaths, names_map, regions_file=None):
    if not all(fp.endswith('.vcf.gz') for fp in fpaths):
        logger.critical('Upset suppots only VCFs')
    sites_file = _run_bcftools_isec(work_dir, fpaths, regions_file)
    upset_file = join(output_dir, 'data.csv')

    with open(sites_file) as f, open(upset_file, 'w') as out:
        out.write('Variant;' + ';'.join(names_map.get(fp, fp) for fp in fpaths) + '\n')
        for l in f:
            chrom, pos, ref, alt, mask = l.split('\t')
            mask = mask.strip()
            mask = map(int, mask)
            out.write(str(chrom) + ':' + str(pos) + ';' + ';'.join(str(b) for b in mask) + '\n')

    return upset_file


def _run_bcftools_isec(work_dir, fpaths, regions_file=None):
    vcfs = ' '.join(fpaths)
    regions_file = ('-T ' + regions_file) if regions_file else ''
    out_dir = join(work_dir, 'bcftools_isec_all')
    cmd = 'bcftools isec {vcfs} -f.,PASS -Oz {regions_file} -p {out_dir}'.format(**locals())
    run(cmd)
    return verify_file(join(out_dir, 'sites.txt'), is_critical=True)


def _intersect_vcfs(work_dir, fpaths, regions_file=None):
    intersection_size_by_subset = dict()
    print('')
    print('Files: ' + str(fpaths))
    print('All subsets:')
    for i in range(1, len(fpaths) + 1):
        for comb in itertools.combinations(fpaths, i):
            print(comb)
            intersection_size_by_subset[tuple(comb)] = 0
    print('')

    sites_file = _run_bcftools_isec(work_dir, fpaths, regions_file)

    with open(sites_file) as f:
        for l in f:
            mask = l.split('\t')[4].strip()
            mask = map(int, mask)
            var_occurences_subset = list(itertools.compress(fpaths, mask))
            for i in range(1, len(var_occurences_subset) + 1):
                for sub_subset in itertools.combinations(var_occurences_subset, i):
                    intersection_size_by_subset[tuple(sub_subset)] += 1

    return intersection_size_by_subset


def _calc_set_intersection(work_dirpath, sorted_files_subset, intersection_file_by_subset):
    """ Intersect using BED tools
    """
    if len(sorted_files_subset) == 1:
        return sorted_files_subset[0]
    if tuple(sorted_files_subset) in intersection_file_by_subset:
        return intersection_file_by_subset[tuple(sorted_files_subset)]

    intersection_file = None
    for fpath in sorted_files_subset:
        remaining_subset = [_f for _f in sorted_files_subset if _f != fpath]  # sorted_beds_subset.index(b) > sorted_beds_subset.index(bed)]
        if not remaining_subset:
            continue
        print('comparing ' + basename(fpath) + ' and ' + str([basename(k) for k in remaining_subset]))
        subset_intersection_file = intersection_file_by_subset.get(tuple(sorted(remaining_subset)))
        if not subset_intersection_file:
            subset_intersection_file = _calc_set_intersection(work_dirpath, remaining_subset, intersection_file_by_subset)
        intersection_file = _intersect_pair(work_dirpath, subset_intersection_file, fpath)
        intersection_file_by_subset[tuple(sorted(remaining_subset))] = subset_intersection_file
    intersection_file_by_subset[tuple(sorted_files_subset)] = intersection_file
    return intersection_file

def _intersect_pair(work_dirpath, file1, file2):
    file1, file2 = sorted([file1, file2])
    output_fpath = join(work_dirpath, splitext(basename(file1))[0] + '__' + basename(file2))
    output_fpath = re.sub(r'.gz$', '', output_fpath)
    if not isfile(output_fpath):
        print('intersect_pair: ' + splitext(basename(file1))[0] + ' and ' + splitext(basename(file2))[0])
        run('bedtools intersect -a ' + file1 + ' -b ' + file2 + ' > ' + output_fpath)
    return output_fpath


def save_venn_diagram_data(size_by_set, names_map):
    data = []
    for venn_set, size in size_by_set.items():
        set_info = dict()
        set_info['size'] = size
        # if isinstance(venn_set, tuple):
        set_info['sets'] = [names_map.get(n, n) for n in venn_set]
        # else:
        #     set_info['sets'] = [venn_set]
        # if isinstance(venn_set, int):
        #     set_info['label'] = label_by_set[venn_set]
        data.append(set_info)
    return json.dumps(sorted(data, key=lambda x: x['sets']))


def write_html(output_dir, json_txt, bed_fpaths):
    output_html = join(output_dir, 'venn.html')

    # def _get_static_file(_fname):
    #     return join(dirname(abspath(reporting.__file__)), 'static', _fname)
    write_static_html_report({
            'title': 'Venn comparison for ' + ', '.join([basename(bf) for bf in bed_fpaths]),
            'diagram_data': json_txt,
        }, output_html,
        tmpl_fpath=join(dirname(abspath(__file__)), 'venn_template.html'),
        extra_js_fpaths=['venn.js', 'd3.min.js', 'd3.tip.js', 'draw_venn_diagram.js'])

    # output_html = join(output_dir, 'venn.html')
    # with open(join(dirname(__file__), 'venn_template.html')) as tmpl_f:
    #     html = tmpl_f.read()
    # html = html.replace('{title}', 'Venn comparison for ' + ', '.join([basename(bf) for bf in bed_fpaths]))
    # html = html.replace('{diagram_data}', json_txt)
    # html = html.replace('{draw_venn_diagram.js}', open(join(dirname(__file__), 'draw_venn_diagram.js')).read())
    # html = html.replace('{venn.js}', open(join(dirname(__file__), 'venn.js')).read())
    # with open(output_html, 'w') as out_f:
    #     out_f.write(html)
    
    return output_html


def check_output(cmdl):
    print(cmdl)
    return subprocess.check_output(cmdl, shell=True)


def bedsize(bed):
    size = check_output("cat " + bed + " | awk -F'\\t' 'BEGIN{ SUM=0 }{ SUM+=$3-$2 }END{ print SUM }'")
    if six.PY3:
        size = size.decode()
    print('Size of ' + basename(bed) + ': ' + size)
    return int(size)


def vcfsize(vcf_file):
    if vcf_file.endswith('.gz'):
        size = check_output("zgrep -v ^# " + vcf_file + " | grep PASS | wc -l")
    else:
        size = check_output("grep -v ^# " + vcf_file + " | grep PASS | wc -l")
    if six.PY3:
        size = size.decode()
    print('Size of ' + basename(vcf_file) + ': ' + size)
    return int(size)

