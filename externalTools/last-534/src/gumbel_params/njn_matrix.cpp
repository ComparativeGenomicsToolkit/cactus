/* $Id: njn_matrix.cpp 183505 2010-02-18 16:10:58Z boratyng $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: njn_matrix.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <ncbi_pch.hpp>

#include "njn_matrix.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(Njn);

   static const MatrixIO::Format FORMAT = MatrixIO::GENERAL;
   static MatrixIO::Format format = MatrixIO::GENERAL;


MatrixIO::Format MatrixIO::getFormat ()
{
   return format;
}

void MatrixIO::setFormat (Format format_)
{
   format = format_;
}

MatrixIO::Format MatrixIO::clearFormat ()
{
   return format = MatrixIO::GENERAL;
}
