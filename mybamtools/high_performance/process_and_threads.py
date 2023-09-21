import time
import os
from multiprocessing import Process
from threading import Thread
from threading import current_thread


def task_p(x: int):
    time.sleep(2)
    print(f"Sub Progress - {os.getpid()}: task {x}: result - {2**x}")


def task_t(x: int):
    time.sleep(2)
    print(f"Sub Thread - {current_thread().name}: task {x}: result - {2**x}")


def running_single_progress():
    print("****** Single Progress ******")
    print(f"Main Progress - {os.getpid()}")
    start = time.time()

    task_p(1)
    task_p(2)

    end= time.time()
    print(f"Main Progress - {os.getpid()}: Use: {end - start:.3f}s")


def running_single_thread():
    print("****** Single Thread ******")
    print(f"Main Thread - {current_thread().name}")
    start = time.time()

    task_t(1)
    task_t(2)

    end= time.time()
    print(f"Main Thread - {current_thread().name}: Use: {end - start:.3f}s")


def running_multi_progress():
    print("****** Multi Progress ******")
    print(f"Main Progress - {os.getpid()}")
    start = time.time()

    p1 = Process(target=task_p, args=(1,)) # tuple
    p2 = Process(target=task_p, args=(2,))
    print("---run sub progress---\n")

    p1.start()
    p2.start()
    p1.join()
    p2.join()

    end = time.time()
    print(f"Main Progress - {os.getpid()}: Use: {end - start:.3f}s")


def running_multi_thread():
    print("****** Multi Thread ******")
    print(f"Main Thread - {current_thread().name}")
    start = time.time()

    t1 = Thread(target=task_t, args=(1,)) # tuple
    t2 = Thread(target=task_t, args=(2,))  # tuple
    print("---run sub thread---\n")

    t1.start()
    t2.start()
    t1.join()
    t2.join()

    end = time.time()
    print(f"Main Thread - {current_thread().name}: Use: {end - start:.3f}s")


if __name__ == "__main__":
    # running_single_progress()
    running_multi_progress()
    # running_single_thread()
    running_multi_thread()

