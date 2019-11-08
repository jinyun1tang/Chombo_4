#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdlib>
#include <iostream>
#include "Chombo_Interval.H"
#include "Chombo_parstream.H"
#include "Chombo_BaseNamespaceHeader.H"

std::ostream&
operator<< (std::ostream&   os,            const Interval& dit)
{
  os << " (" << dit.m_begin << "," << dit.m_end << ") " ;
  return os;
}
#include "Chombo_BaseNamespaceFooter.H"
