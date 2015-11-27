/*
 * LatticeMatch calculator - class used for angles
 * Copyright (C) 2015 Andreas Grois
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
    fullcircle= false;
    sortby = SRT_SIZE;
}

anglerange::anglerange(const angleclass &lower, const angleclass &upper)
{
    lowerborder = lower;
    upperborder = upper;
    upperset = true;
    lowerset = true;
    fullcircle = false;
    sortby = SRT_SIZE;
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
    fullcircle = false;
//    lowerborder = 0.0;
//    upperborder = 0.0;
//commented out, as users (I) should not rely on this.
}

void anglerange::setsorttype(const anglerange::sorttype value)
{
    sortby=value;
}

const anglerange::sorttype anglerange::getsorttype()
{
    return(sortby);
}

bool anglerange::iscircle() const
{
    return(fullcircle && !isempty());
}

void anglerange::setcircle(bool value)
{
    if(value)
    {
        lowerset=true;
        upperset=true;
        lowerborder=0.0;
        upperborder=0.0;
    }
    fullcircle=value;
}

bool anglerange::isinside(const angleclass &val) const
{
    if(isempty())
    {
        return(false);
    }
    else
    {
        if(fullcircle)
        {
            return(true);
        }
        else if(lowerborder>upperborder)
        {
            return((val>=lowerborder || val<=upperborder));
        }
        else
        {
            return(val>=lowerborder && val<=upperborder);
        }
    }
}

anglerange anglerange::overlap(const anglerange &other) const
{
    //storage for the return value:
    anglerange retval; //anglerange constructor without arguments: emtpy range [0:0]
    retval.setsorttype(sortby); //return value inherits current sorttype.
    //first: if one of the input ranges is empty, don't do anything, leave retval empty!

    if(!(isempty()) && !(other.isempty()))
    {
        //special treatment if one of the ranges is a full circle:
        if(fullcircle)
        {
            retval = other;
            retval.setsorttype(sortby);
        }
        else if(other.fullcircle)
        {
            retval = *this;
        }
        //Here we need to take care to always follow the counter-clockwise convention.
        //This leaves us with four possibilities:
        //a) Both ranges contain the zero-line
        //b) This range contains the zero-line, the other doesn't
        //c) Vice versa
        //d) Both ranges are regular.

        //first check if this range contains the zero-line
        else if(lowerborder>upperborder)
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

anglerange anglerange::combine(const anglerange &other) const
{
    anglerange retval;
    retval.setsorttype(sortby); //return value inherits current sorttype.
    //check if any of the ranges is empty. If yes: return an empty range as well.
    if(!(isempty()) && !(other.isempty())){
        //again: special treatment if one of the ranges is a full circle:
        if(fullcircle)
        {
            retval = *this;
        }
        if(other.fullcircle)
        {
            retval = other;
            retval.setsorttype(sortby);
        }
        //Here we need to take care to always follow the counter-clockwise convention.
        //This leaves us with four possibilities:
        //a) Both ranges contain the zero-line
        //b) This range contains the zero-line, the other doesn't
        //c) Vice versa
        //d) Both ranges are regular.

        //first check if this range contains the zero-line
        else if(lowerborder>upperborder)
        {
            //ok, this range is around zero. Let's check the other one too:
            if(other.getlower()>other.getupper())
            {
                //ok, both ranges contain the zero-line.
                //let's see if the result is a full circle
                if(lowerborder<=other.upperborder || other.lowerborder <=upperborder)
                {
                    //full circle
                    retval.setcircle(true);
                }
                else
                {
                    retval.setlower(fmin(lowerborder.getval(),other.getlower().getval()));
                    retval.setupper(fmax(upperborder.getval(),other.getupper().getval()));
                }
            }
            else
            {
                //this contains zero, other doesn't.
                //This cannot fully lie within other, but other can lie fully within this.
                //The result does obviously contain zero.
                //Let's check if other is within this.
                if(other.getupper()<=upperborder || other.getlower()>=lowerborder){
                    retval.setlower(lowerborder);
                    retval.setupper(upperborder);
                }
                //again: check if result is a full circle
                else if(other.lowerborder<=upperborder && other.upperborder>=lowerborder)
                {
                    retval.setcircle(true);
                }
                else if(other.getupper()>=lowerborder) //not >, as per definition upperborder is inside range
                {
                    retval.setlower(other.getlower());
                    retval.setupper(upperborder);
                }
                else if(other.getlower()<=upperborder) //not <, as again upperborder is inside.
                {
                    retval.setlower(lowerborder);
                    retval.setupper(other.getupper());
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
                    retval.setlower(other.getlower());
                    retval.setupper(other.getupper());
                }
                //full circle?
                else if(lowerborder <= other.getupper() && upperborder >= other.getlower())
                {
                    retval.setcircle(true);
                }
                else if(lowerborder <= other.getupper())
                {
                    retval.setlower(other.getlower());
                    retval.setupper(upperborder);
                }
                else if(upperborder >= other.getlower())
                {
                    retval.setlower(lowerborder);
                    retval.setupper(other.getupper());
                }
            }
            else //both ranges regular
            {
                //this also means: none of them contains zero, the result cannot be a full circle
                //check if there's an overlap
                double curmax = fmin(other.getupper().getval(),upperborder.getval());
                double curmin = fmax(other.getlower().getval(),lowerborder.getval());
                if(curmax>=curmin){
                    //there is an overlap - set limits
                    retval.setlower(fmin(other.getlower().getval(),lowerborder.getval()));
                    retval.setupper(fmax(other.getupper().getval(),upperborder.getval()));
                }
            }
        }
    }
    return(retval);
}

bool anglerange::operator<(const anglerange &other) const
{
    if(isempty() || other.isempty())
    {
        return(false);
    }
    else
    {
        switch(sortby)
        {
        case SRT_LOWER:
            return(lowerborder < other.lowerborder);
            break;
        case SRT_UPPER:
            return(upperborder < other.upperborder);
            break;
        case SRT_SIZE:
            //I know that this if elseif thing can be written as a single statement. This is for readability.
            if(other.fullcircle)
            {
                return(!fullcircle);
            }
            else if(fullcircle)
            {
                return(false);
            }
            else
            {
                //angleclass takes care of zero crossing. C++ is awesome.
                return(upperborder - lowerborder < other.upperborder - other.lowerborder);
            }
            break;
        default:
            throw std::out_of_range("Invalid sort field set for anglerange! Go, debug.\n");
        }
    }
}

bool anglerange::operator>(const anglerange &other) const
{
    //as == doesn't care about sort order but fully compares, we have to implement two sort operators fully...
    if(isempty() || other.isempty())
    {
        return(false);
    }
    else
    {
        switch(sortby)
        {
        case SRT_LOWER:
            return(lowerborder > other.lowerborder);
            break;
        case SRT_UPPER:
            return(upperborder > other.upperborder);
            break;
        case SRT_SIZE:
            if(fullcircle)
            {
                return(!other.fullcircle);
            }
            else if(other.fullcircle)
            {
                return(false);
            }
            else
            {
                //angleclass takes care of zero crossing. C++ is awesome.
                return(upperborder - lowerborder > other.upperborder - other.lowerborder);
            }
            break;
        default:
            throw std::out_of_range("Invalid sort field set for anglerange! Go, debug.\n");
        }
    }
}

bool anglerange::operator<=(const anglerange &other) const
{
    if(isempty() || other.isempty())
    {
        return(false);
    }
    else
    {
        return(!(operator>(other)));
    }
}
bool anglerange::operator>=(const anglerange &other) const
{
    if(isempty() || other.isempty())
    {
        return(false);
    }
    else
    {
        return(!(operator<(other)));
    }
}


bool anglerange::operator==(const anglerange &other) const
{
    if(!(isempty()) && !(other.isempty()))
    {
        //both ranges set, we need to compare field by field
        return( (lowerborder == other.lowerborder) &&  (upperborder == other.upperborder) );
    }
    else
        return(isempty() && other.isempty());
}
bool anglerange::operator!=(const anglerange &other) const
{
    //needn't check for isempty here, as == already checks for it.
    return(!operator==(other));
}
