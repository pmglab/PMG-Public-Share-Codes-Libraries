# 批量检查指定GDC manifest文件中对应测序数据文件的md5值是否匹配

import os
import pandas as pd
import argparse
from multiprocessing import Pool

# 指定输入参数与信息
parser=argparse.ArgumentParser(description="Batch check md5 for given sequence file within manifest.")
parser.add_argument('--manifest','-m',help='Location of manifest file.')
parser.add_argument('--check_folder', '-f',help='Location of folder containing GDC sequence data.')
parser.add_argument('--num_threads', '-nt',help='Number of threads concurrently used to check md5 value.')
parser.add_argument("--output_log", '-o', help="Log file of checking results")

args=parser.parse_args()

# 获取manifest信息
manifest_info = pd.read_table(args.manifest)

def md5_status_getter(file_loc, given_md5_value, file_object):
    real_md5 = os.popen(f'md5sum {file_loc}').read()
    if real_md5.split()[0] == given_md5_value:
        file_object.write(f"{file_loc} had been verified successfully.\n")
    else:
        file_object.write(f"\nFail warning: Be advised, {file_loc} failed md5 check.\n")


with open(args.output_log, 'w') as f:
    print('Parent process %s.' % os.getpid())
    p = Pool(int(args.num_threads))
    for i in manifest_info.index:
        file_loc = os.path.join(args.check_folder, manifest_info['id'][i], manifest_info['filename'][i])
        given_md5_value = manifest_info['md5'][i]
        p.apply_async(md5_status_getter, args=(file_loc, given_md5_value, f,))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')