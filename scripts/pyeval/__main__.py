from pprint import pprint

from args import parse_args

# ['--executable', 'camel',
#          '--args', 
#          '"--threads 32 --comment testing"',
#          '--database', './sample.db',
#          '--work-dir', '/tmp/camel']

try:
    args = parse_args()
    pprint(args)
except Exception as e:
    print(e)
