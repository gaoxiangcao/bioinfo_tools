#!/usr/bin/env python
# _*_ coding: UTF-8 _*_
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-a', '--input', help='x', type=str, required=True,default='stdout')


ARGS = parser.parse_args()

print(ARGS.input)



