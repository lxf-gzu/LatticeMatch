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

#ifndef ANGLERANGE_H
#define ANGLERANGE_H

#include "angleclass.h"

class anglerange
{
private:
    angleclass lowerborder;
    angleclass upperborder;
    bool upperset;
    bool lowerset;
public:
    //Default constructor: Both limits 0.0, marked as empty
    anglerange();
    anglerange(const angleclass &lower, const angleclass &upper);
    void setlower(const angleclass &newlower);
    void setupper(const angleclass &newupper);
    void setempty(); //marks the range as empty.

    bool isempty() const; //returns true when the range is empty (or when one limit isn't set)
    angleclass getlower() const; //warning: Does not check if empty
    angleclass getupper() const; //warning: Does not check if empty

    //this function calculates the overlap between this anglerange and another one, returning it as
    //a new anglerange.
    anglerange overlap(const anglerange &other);

};

#endif // ANGLERANGE_H
