# Make.sh = update Makefile.lib or Makefile.list or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.list

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

# create individual style files
# called by "make machine"

if (test $1 = "Makefile.lib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.lib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.lib

# edit Makefile.shlib, for creating shared lib
# called by "make makeshlib"
# use current list of *.cpp and *.h files in src dir w/out main.cpp

elif (test $1 = "Makefile.shlib") then

  list=`ls -1 *.cpp | sed s/^main\.cpp// | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.shlib
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.shlib

# edit Makefile.list
# called by "make makelist"
# use current list of *.cpp and *.h files in src dir

elif (test $1 = "Makefile.list") then

  list=`ls -1 *.cpp | tr "[:cntrl:]" " "`
  sed -i -e "s/SRC =	.*/SRC =	$list/" Makefile.list
  list=`ls -1 *.h | tr "[:cntrl:]" " "`
  sed -i -e "s/INC =	.*/INC =	$list/" Makefile.list

fi
