#!/usr/bin/env python 
import os
import subprocess

###Retrieve version number from git describe
### Copied with modifications from: http://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
def git_version():
 def _minimal_ext_cmd(cmd):
  # construct minimal environment
  env = {}
  for k in ['SYSTEMROOT', 'PATH']:
   v = os.environ.get(k)
   if v is not None:
    env[k] = v
  out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
  return out

 try:
  out = _minimal_ext_cmd(['git', 'describe'])
  GIT_REVISION = out.strip().decode('ascii').split("-")[0]
 except OSError:
  GIT_REVISION = "Unknown"
 return GIT_REVISION
