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
 * This is a small class that doesn't do much more than keep its contents in the range between
 * [0:2*pi[. It is a part of the LatticeMatch program.
 *
 * To contact the author either use electronic mail: Andreas Grois <andreas.grois@jku.at>
 * or write to:
 * Andreas Grois, Institute for Semiconductor and SolidState physics, Johannes Kepler University,
 * Altenbergerstraße 69, 4040 Linz, AUSTRIA
 */

#include "angleclass.h"
#include <cmath>


double angleclass::shiftinrange(const double number)
{
    //Shifts number back in the range between zero and 2 pi
    return(number-2.0*M_PI*floor(number/(2.0*M_PI)));
}

angleclass::angleclass()
{
    //Default constructor. Sets value to zero - what else?
    value=0.0;
}

angleclass::angleclass(const double setme)
{
    value=shiftinrange(setme);
}

void angleclass::setval(const double setme)
{
    value=shiftinrange(setme);
}

double angleclass::getval() const
{
    return(value);
}

angleclass angleclass::operator *(const angleclass other) const
{
    return(angleclass(value*other.value));
}


angleclass angleclass::operator-(const angleclass other) const
{
    return(angleclass(value-other.value));
}

angleclass angleclass::operator+(const angleclass other) const
{
    return(angleclass(value+other.value));
}

angleclass angleclass::operator/(const angleclass other) const
{
    return(angleclass(value/other.value));
}

bool angleclass::operator<(const angleclass other) const
{
    return(value<other.value);
}

bool angleclass::operator>(const angleclass other) const
{
    return(value>other.value);
}

bool angleclass::operator<=(const angleclass other) const
{
    return(value<=other.value);
}

bool angleclass::operator>=(const angleclass other) const
{
    return(value>=other.value);
}

bool angleclass::operator==(const angleclass other) const
{
    return(value==other.value);
}

bool angleclass::operator!=(const angleclass other) const
{
    return(value!=other.value);
}
