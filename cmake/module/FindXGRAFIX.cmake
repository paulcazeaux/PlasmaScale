# - Find XGRAFIX
# Find the native XGRAFIX includes and library
#
#  XGRAFIX_INCLUDES    - where to find xgrafix.h
#  XGRAFIX_LIBRARIES   - List of libraries when using XGRAFIX.
#  XGRAFIX_FOUND       - True if XGRAFIX found.

if (XGRAFIX_INCLUDES)
  # Already in cache, be silent
  set (XGRAFIX_FIND_QUIETLY TRUE)
endif (XGRAFIX_INCLUDES)

find_path (XGRAFIX_INCLUDE xgrafix.h)

find_library (XGRAFIX_LIBRARY NAMES XGC250)
find_library(TCL_LIBRARY NAMES tcl8.6 HINTS "/usr/local/opt/tcl-tk/lib")
find_library(TK_LIBRARY NAMES tk8.6 HINTS "/usr/local/opt/tcl-tk/lib")
find_library(XPM_LIBRARY NAMES Xpm)
find_library(X11_LIBRARY NAMES X11)
find_library(M_LIBRARY NAMES m)
find_library(DL_LIBRARY NAMES dl)

set(XGRAFIX_INCLUDES ${XGRAFIX_INCLUDE} "/opt/X11/lib")
set(XGRAFIX_LIBRARIES ${XGRAFIX_LIBRARY} ${TCL_LIBRARY} ${TK_LIBRARY} ${XPM_LIBRARY} ${X11_LIBRARY} ${M_LIBRARY} ${DL_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set XGRAFIX_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (XGRAFIX DEFAULT_MSG XGRAFIX_LIBRARIES XGRAFIX_INCLUDES)

mark_as_advanced (XGRAFIX_LIBRARIES XGRAFIX_INCLUDES)