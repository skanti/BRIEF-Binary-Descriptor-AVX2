divert(`-1')
# forloop(var, from, to, body) - simple version
define(`forloop', `pushdef(`$1', `$2')_forloop($@)popdef(`$1')')
define(`_forloop', `$4`'ifelse($1, `$3', `', `define(`$1', incr($1))$0($@)')')

# m4_var(var, number) - gets numbered variable
define(`m4_var',$1_$2)

divert`'dnl

