params.root_folder = "/media/davidp/Data/coding/nf-testing/nf-testing"

process foo {
  debug true
  script:
  """
  echo foo task path: \$PWD
  """ 
}

process bar {
  debug true
  script:
  '''
  echo bar task path: $PWD
  ''' 
}

workflow {
  foo()
  bar()
}