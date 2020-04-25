#!/usr/bin/python3

#def test_sys_module():
#  mod = 'sys'
#  try:
#    import sys
#    print("{:>15} Ok".format(mod))
#    return(1)
#  except:
#    print("{:>15} Failed".format(mod))
#    return(0)

def test_importlib_module():
  mod = 'importlib'
  try:
    import importlib
    print("{:>15} Ok".format(mod))
    return(1)
  except:
    print("{:>15} Failed".format(mod))
    sys.exit(-1)
    return(0)

def test_any_python_module(mod):
  try:
    test_import = importlib.import_module(mod)
    print("{:>15} Ok".format(mod))
    return(1)
  except:
    print("{:>15} Failed".format(mod))
    return(0)

print("Test all modules required for \'md2nmr\'")
status = 1
if test_importlib_module():
  import importlib
  status *= test_any_python_module('numpy')
  status *= test_any_python_module('prettytable')
  status *= test_any_python_module('argparse')
  status *= test_any_python_module('subprocess')
  status *= test_any_python_module('matplotlib')
  status *= test_any_python_module('lmfit')
  status *= test_any_python_module('tqdm')
else:
  print("Tets could not be proceeded. Please install \'importlib\' python module and start again")
  status = 0

print("Overall status = {}".format(status))

