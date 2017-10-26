#!/usr/bin/env python
from __future__ import print_function

import re
import subprocess
import json

import six
from ngs_utils.reporting.reporting import write_static_html_report
from os.path import dirname, basename, join, splitext, isfile, abspath, pardir


def call(cmdl):
    print(cmdl)
    subprocess.call(cmdl, shell=True)


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


total_calls = 0


def intersect_pair(work_dirpath, file1, file2):
    file1, file2 = sorted([file1, file2])
    output_fpath = join(work_dirpath, splitext(basename(file1))[0] + '__' + basename(file2))
    output_fpath = re.sub(r'.gz$', '', output_fpath)
    if not isfile(output_fpath):
        print('intersect_pair: ' + splitext(basename(file1))[0] + ' and ' + splitext(basename(file2))[0])
        call('bedtools intersect -a ' + file1 + ' -b ' + file2 + ' > ' + output_fpath)
    return output_fpath


def calc_set_intersection(work_dirpath, sorted_files_subset, intersection_file_by_subset):
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
            subset_intersection_file = calc_set_intersection(work_dirpath, remaining_subset, intersection_file_by_subset)
        intersection_file = intersect_pair(work_dirpath, subset_intersection_file, fpath)
        intersection_file_by_subset[tuple(sorted(remaining_subset))] = subset_intersection_file
    intersection_file_by_subset[tuple(sorted_files_subset)] = intersection_file
    return intersection_file


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


def run(work_dirpath, fpaths):
    intersection_file_by_subset = dict()
    calc_set_intersection(work_dirpath, sorted(fpaths), intersection_file_by_subset)

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

