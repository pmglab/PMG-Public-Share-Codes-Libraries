# 批量检查指定GDC manifest文件中对应测序数据文件的md5值是否匹配

import os
import pandas as pd
import argparse
from multiprocessing import Pool
from multiprocessing import Manager

# 指定输入参数与信息
parser=argparse.ArgumentParser(description="Batch check md5 for given sequence file within manifest.")
parser.add_argument('--manifest','-m',help='Location of manifest file.')
parser.add_argument('--check_folder', '-f',help='Location of folder containing GDC sequence data.')
parser.add_argument('--num_threads', '-nt',help='Number of threads concurrently used to check md5 value. (recommend <=10, IO sensitive task)')
parser.add_argument("--output_log", '-o', help="Log file of checking results")

args=parser.parse_args()

# 获取manifest信息
manifest_info = pd.read_table(args.manifest)

def md5_status_getter(file_loc, given_md5_value, return_list):
    print(f"Start to check {file_loc} md5 status.")
    real_md5 = os.popen(f'md5sum {file_loc}').read()
    if real_md5.split()[0] == given_md5_value:
        return_list.append(f"{file_loc} had been verified successfully.")
    else:
        return_list.append(f"\nFail warning: Be advised, {file_loc} failed md5 check. Provided md5 was {given_md5_value}, whereas calculated md5 was {real_md5.split()[0]}\n")

if __name__=='__main__':
    # Create a list to store md5-check results
    manager = Manager()
    return_list = manager.list()
    # Start multi-processing
    print('Parent process %s.' % os.getpid())
    p = Pool(int(args.num_threads))
    for i in manifest_info.index:
        file_loc = os.path.join(args.check_folder, manifest_info['id'][i], manifest_info['filename'][i])
        given_md5_value = manifest_info['md5'][i]
        p.apply_async(md5_status_getter, args=(file_loc, given_md5_value, return_list))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')
    # Write result list into log file
    with open(args.output_log, 'w') as f: 
        f.write("\n".join(list(return_list)))