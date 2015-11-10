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

#include "angleset.h"

angleset::angleset()
{
    //storage should be initialized as an empty vector, all we need here is:
    consistent=true;
}

angleset::angleset(const anglerange &firstrange)
{
    assert(storage.size()==0);
    if(!firstrange.isempty())
    {
        storage.push_back(firstrange);
        storage[0].setsorttype(anglerange::SRT_LOWER);
    }
    consistent=true;
}

angleset::angleset(const angleclass &firstlower, const angleclass &firstupper)
{
    assert(storage.size()==0);
    storage.push_back(anglerange(firstlower,firstupper));
    storage[0].setsorttype(anglerange::SRT_LOWER);
    consistent=true;
}

void angleset::reserve(size_t n)
{
    storage.reserve(n);
}

void angleset::combine()
{
    //the following should hopefully also work when storage is empty.
    if(storage.size()>=2){ //only do something if there are at least two elements
        for(std::vector<anglerange>::iterator curpos=storage.end()-2; //second last element
            curpos>=storage.begin(); //learned the hard way, that for an empty vector .begin() returns (unsigned)0, so .begin()-2>.begin()...
            curpos--){
            std::vector<anglerange>::iterator compare = storage.end()-1;
            while(compare!=curpos)
            {
                anglerange tmp = curpos->combine(*compare);
                if(!(tmp.isempty()))
                {
                    *curpos=tmp;
                    compare--;  //I'm not certain if it's necessary to first decrement before deleting,
                                //but better safe than sorry.
                    storage.erase(compare+1);
                }
                else
                    compare--;
            }
        }
    }
    consistent=true;
}

bool angleset::isempty() const
{
    return(storage.empty());
}

bool angleset::iscircle()
{
    if(!consistent)
        combine();
    //if only one element in storage should be a full circle, combine will reduce storage to a single element vector.
    if(!storage.empty())
    {
        return(storage.front().iscircle());
    }
    else{
        return(false);
    }
}

void angleset::add(const anglerange &value)
{
    //quick and dirty: We add the range, sort through, and combine if necessarry.
    //this means: we *should not reuse* this function for adding anglesets.
    if(!value.isempty())
    {
        storage.push_back(value);
        storage.back().setsorttype(anglerange::SRT_LOWER);
        combine();
    }
}

void angleset::add(const angleset &value)
{
    storage.reserve(storage.size()+value.storage.size());
    storage.insert(storage.end(),value.storage.begin(),value.storage.end());
    //if value is a valid angleset, sort order is set to SRT_LOWER for all fields.
    combine();
}

void angleset::add(const double &lower, const double &upper)
{
    storage.push_back(anglerange(lower,upper));
    storage.back().setsorttype(anglerange::SRT_LOWER);
    combine();
}


angleset angleset::overlap(const anglerange &other)
{
    if(!consistent)
        combine();
    angleset retval;
    for_each(
                storage.begin(),
                storage.end(),
                [&](anglerange i)
                {
                    retval.add(i.overlap(other)); //given that add checks if the value passed to it is empty: this should work?!?
                }
    );
    return(retval);
}

angleset angleset::overlap(const angleset &other)
{
    //just add the overlap of each anglerange in other ;-)
    angleset retval;
    for_each(
                other.storage.begin(),
                other.storage.end(),
                [&](anglerange i)
                {
                    retval.add(overlap(i));
                }
                );
    return(retval);
}

std::vector<anglerange> angleset::getranges()
{
    return storage;
}

const std::vector<anglerange>& angleset::getrangesref()
{
    return storage;
}

void angleset::clear()
{
    storage.clear();
}

void angleset::sort()
{
    std::sort(storage.begin(),storage.end());
}
