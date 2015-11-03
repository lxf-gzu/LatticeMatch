/*
 * LatticeMatch calculator - class used for angles
 * Copyright (C) 2012 Andreas Grois
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 *
 * A small class that deals with ranges of angles. Its members are two angles, maximum and minimum.
 * The main purpose of this class is to provide an easy way to determine the overlap of two ranges.
 * There is no support for disjoint ranges (yet).
 * Probably that would better be left for a separate class anyhow. Think: vector of anglerange objects
 *
 * There is an important convention with this class:
 * Ranges are always clockwise from the lower to the upper bound.
 * They include the lower and the upper bound, as this is required to be able to contain individual points.
 *
 * This class is part of the LatticeMatch program.
 *
 * To contact the author either use electronic mail: Andreas Grois <andreas.grois@jku.at>
 * or write to:
 * Andreas Grois, Institute for Semiconductor and SolidState physics, Johannes Kepler University,
 * Altenbergerstraße 69, 4040 Linz, AUSTRIA
 */

#include "anglerange.h"
#include <cmath>

anglerange::anglerange()
{
    lowerborder = angleclass(0.0);
    upperborder = angleclass(0.0);
    upperset = false;
    lowerset = false;
}

anglerange::anglerange(const angleclass &lower, const angleclass &upper)
{
    lowerborder = lower;
    upperborder = upper;
    upperset = true;
    lowerset = true;
}

bool anglerange::isempty() const
{
    return(!upperset || !lowerset);
}

angleclass anglerange::getupper() const
{
    return(upperborder);
}

angleclass anglerange::getlower() const
{
    return(lowerborder);
}

void anglerange::setupper(const angleclass &newupper)
{
    upperborder=newupper;
    upperset = true;
}

void anglerange::setlower(const angleclass &newlower)
{
    lowerborder=newlower;
    lowerset = true;
}

void anglerange::setempty()
{
    upperset = false;
    lowerset = false;
    lowerborder = 0.0;
    upperborder = 0.0;
}

anglerange anglerange::overlap(const anglerange &other)
{
    //storage for the return value:
    anglerange retval; //anglerange constructor without arguments: emtpy range [0:0]
    //first: if one of the input ranges is empty, don't do anything, leave retval empty!
    if(!(isempty()) && !(other.isempty()))
    {
        //Here we need to take care to always follow the counter-clockwise convention.
        //This leaves us with four possibilities:
        //a) Both ranges contain the zero-line
        //b) This range contains the zero-line, the other doesn'T
        //c) Vice versa
        //d) Both ranges are regular.

        //first check if this range contains the zero-line
        if(lowerborder>upperborder)
        {
            //ok, this range is around zero. Let's check the other one too:
            if(other.getlower()>other.getupper())
            {
                //ok, both ranges contain the zero-line. This makes some things easier, as we know there's
                //an overlap
                retval.setlower(fmax(lowerborder.getval(),other.getlower().getval()));
                retval.setupper(fmin(upperborder.getval(),other.getupper().getval()));
            }
            else
            {
                //this contains zero, other doesn't.
                //This cannot fully lie within other, but other can lie fully within this.
                //The result does obviously not contain zero.
                //Let's check if other is within this.
                if(other.getupper()<=upperborder || other.getlower()>=lowerborder){
                    retval.setlower(other.getlower());
                    retval.setupper(other.getupper());
                }
                else if(other.getupper()>=lowerborder) //not >, as per definition upperborder is inside range
                {
                    retval.setlower(lowerborder);
                    retval.setupper(other.getupper());
                }
                else if(other.getlower()<=upperborder) //not <, as again upperborder is inside.
                {
                    retval.setlower(other.getlower());
                    retval.setupper(upperborder);
                }
            }
        }
        else
        {
            //this is regular, is other as well?
            if(other.getlower()>other.getupper())
            {
                //so, other cannot lie in this, but this can be in other.
                //Let's check for that.
                if(upperborder<=other.getupper() || lowerborder>=other.getlower()){
                    retval.setlower(lowerborder);
                    retval.setupper(upperborder);
                }
                else if(lowerborder <= other.getupper())
                {
                    retval.setlower(lowerborder);
                    retval.setupper(other.getupper());
                }
                else if(upperborder >= other.getlower())
                {
                    retval.setlower(other.getlower());
                    retval.setupper(upperborder);
                }
            }
            else //both ranges regular
            {
                double curmax = fmin(other.getupper().getval(),upperborder.getval());
                double curmin = fmax(other.getlower().getval(),lowerborder.getval());
                if(curmax>=curmin){
                    retval.setlower(curmin);
                    retval.setupper(curmax);
                }
            }
        }
    }
    return(retval);
}
