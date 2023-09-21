"""Pool class
apply & apply_async
map & map_async
close, join, terminate
"""

import time
import os
from multiprocessing import Pool


def task_p(x: int):
    time.sleep(2)
    print(f"Sub Progress - {os.getpid()}: task {x}: Result - {2**x}")


def running_multi_progress():
    print("****** Multi Progress ******")
    print(f"Main Progress - {os.getpid()}")
    start = time.time()

    pool = Pool(processes=os.cpu_count() - 4)
    # for x in range(1, 101, 1):
    #     # pool.apply(func=task_p, args=(x,))
    #     pool.apply_async(func=task_p, args=(x,))
    # pool.map(func=task_p, iterable=range(1,3,1))
    pool.map_async(func=task_p, iterable=range(1, 3, 1))
    print("---run sub progress---\n")


    pool.close()
    pool.join()

    end = time.time()
    print(f"Main Progress - {os.getpid()}: Use: {end - start:.3f}s")


if __name__ == "__main__":
    running_multi_progress()

