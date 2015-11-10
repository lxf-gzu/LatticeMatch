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
 * This class extends anglerange with the ability to deal with disjoint ranges.
 *
 * The internal storage isn't exposed, and there is no way to address an individual sub-range.
 * The reason for this is that the functions to add or remove ranges will resort the contents
 * of the internal storage arbitrarily.
 *
 * This class is part of the LatticeMatch program.
 *
 * To contact the author either use electronic mail: Andreas Grois <andreas.grois@jku.at>
 * or write to:
 * Andreas Grois, Institute for Semiconductor and SolidState physics, Johannes Kepler University,
 * Altenbergerstraße 69, 4040 Linz, AUSTRIA
 */

#ifndef ANGLESET_H
#define ANGLESET_H
#include <vector>
#include <cassert>
#include <algorithm>
#include "anglerange.h"

class angleset
{
private:
    std::vector<anglerange> storage;
    bool consistent;
    void combine();
public:
    angleset();
    angleset(const anglerange &firstrange);
    angleset(const angleclass &firstlower, const angleclass &firstupper);

    //returns if the set is empty. Currently this means: internal storage is an empty vector.
    bool isempty() const;
    bool iscircle();//returns true when the whole angleset is a circle.

    //add or remove single ranges
    void add(const anglerange &value);
    //void remove(const anglerange &value);  //not implemented yet
    void add(const double &lower, const double &upper);
    //void remove(const double &lower, const double &upper); //not implemented yet

    //here it gets interesting: add or remove complete sets.
    void add(const angleset &value);
    //void remove(const angleset &value); //not implemented yet

    //in planning: add_deferred functions, that set consistent to false, and do not call combine() at the end.
    //but for this of course each and every other function has to check for consistency

    //having this return a new angleset is not consistent with the above add/remove functions, but it
    //is consistent to older code in anglerange.h
    angleset overlap(const anglerange &other);
    angleset overlap(const angleset &other);

    void reserve(size_t n); //just forwards storage's reserve function.

    //empties the range
    void clear();

    //sorts the range
    void sort();

    //*generates* a new vector consisting of non-overlapping, unique angleranges.
    std::vector<anglerange> getranges();
    //gives you a reference to the internal storage, but const
    //reference, to keep this relatively fast, const so you don't accidentally modify it.
    //If the internal storage format changes, this might be changed to give a copy instead.
    const std::vector<anglerange>& getrangesref(); //I hope this does what I think it does...

};

#endif // ANGLESET_H
