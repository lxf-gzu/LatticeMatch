/*
 * LatticeMatch calculator
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
 * This is a small script that calculates possible lattice matches as discussed in:
 *
 * Fritz, Torsten: Molecular Architecture in Heteroepitaxially Grown Organic Thin Films
 * Dresden: sfps - Wissenschaftlicher Fachverlag, 1999
 * ISBN 3-934264-50-6
 *
 * The advantage of this script compared to the original work of T. Fritz is, that this
 * program does not use a brute-force approach, but instead determines ranges of possible
 * values for the angle between the first lattice vectors of the interface unit cell of
 * substrate and adlayer analytically.
 * To do this, the user has to enter the substrate interface unit cell, giving the length
 * of the lattice vectors a1 and a2, and also the angle between them, alpha.
 * For the interface unit cell of the adlayer the user has to input ranges in which the
 * magnitude of the lattice vectors b1 and b2 may lie, as well as a range for the angle
 * between them, beta. The program then uses these inequalities together with the
 * underdefined set of equations describing a coincident lattice match to determine ranges
 * in which the angle theta between a1 and b1 may lie.
 *
 * According to the work quoted above, a coincident lattice match occurs when both elements
 * in one of the columns of the epitaxy matrix are integer, while in case all entries are
 * integer, one is dealing with a commensurate lattice match.
 *
 * The epitaxy matrix, as given in the work quoted above, reads:
 *
 * ( px, qy )
 * ( qx, py )
 *
 * with:
 * px=b1*sin(alpha-theta)/(a1*sin(alpha))
 * qx=b2*sin(alpha-theta-beta)/(a1*sin(alpha))
 * qy=b1*sin(theta)/(a2*sin(alpha))
 * py=b2*sin(theta+beta)/(a2*sin(alpha))
 *
 * There are two reasons for the length of this program:
 * o) asin isn't unique.
 * o) ranges of angles are a pain to deal with.
 * Both points together made writing the correct conditions for upper and lower bounds
 * of the result ranges a pain, and I'm pretty certain there might still be the one or the other
 * bug hidden.
 *
 * To contact the author either use electronic mail: Andreas Grois <andreas.grois@jku.at>
 * or write to:
 * Andreas Grois, Institute for Semiconductor and SolidState physics, Johannes Kepler University,
 * Altenbergerstraße 69, 4040 Linz, AUSTRIA
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include "anglerange.h"

using namespace std;

enum commargnames{A1,A2,ALPHA,B1MIN,B1MAX,B2MIN,B2MAX,BETAMIN,BETAMAX};

int main(int argc, char* argv[])
{
    //This is largely a copy of what I already did in plain C.
    //There will be a lot of plain C code here, so don't look too close...
    if(argc!=10)
    {
        cout << "Usage: " << argv[0] << " a1 a2 alpha b1min b1max b2min b2max betamin betamax" << std::endl << "Please input angles in degrees." << std::endl;
        return(-1);
    }
    else
    {
        //read in command line arguments
        double commargs[9], swapper;
        unsigned int i,j;
        for(i=1;i<(unsigned)argc;i++)
        {
            //what is this stringstreams stuff everyone is hyped about?!?
           sscanf(argv[i],"%lf",&(commargs[i-1]));
        }
        //Never trust the user. I know that I'll probably be the only one to use this program, but that's just another reason...
        //to make sure that the minimum values are actually smaller than the maximum numbers. Also, we need radians, not degrees.
        if(commargs[B1MIN]*commargs[B2MIN]<0 || commargs[B1MAX]*commargs[B2MAX]<0 || commargs[B1MIN]*commargs[B1MAX]<0)
        {
            std::cerr << "Warning: negative values for b1, b2 don't make any sense. Putting them back in order." << std::endl;
            commargs[BETAMIN]=180-commargs[BETAMIN];
            commargs[BETAMAX]=180-commargs[BETAMAX];
        }
        if(commargs[A1]*commargs[A2]<0)
        {
            std::cerr << "Warning: negative values for a1, a2 don't make any sense. Putting them back in order." << std::endl;
            commargs[ALPHA]=180-commargs[ALPHA];
        }
        for(i=B1MIN;i<=B2MAX;i++)
        {
            commargs[i]=fabs(commargs[i]);
        }
        commargs[A1]=fabs(commargs[A1]);
        commargs[A2]=fabs(commargs[A2]);

        commargs[BETAMIN]=(commargs[BETAMIN]-360*floor(commargs[BETAMIN]/360.0))*M_PI/180.0;
        commargs[BETAMAX]=(commargs[BETAMAX]-360*floor(commargs[BETAMAX]/360.0))*M_PI/180.0;
        commargs[ALPHA]=(commargs[ALPHA]-360*floor(commargs[ALPHA]/360.0))*M_PI/180.0;

        for(i=0;i<3;i++)
        {
            if(commargs[2*i+3]>commargs[2*i+1+3])
            {
                swapper = commargs[2*i+3];
                commargs[2*i+3] = commargs[2*i+1+3];
                commargs[2*i+1+3] = swapper;
            }
        }
        if(commargs[BETAMAX]-commargs[BETAMIN]>M_PI)
        {
            std::cerr << "Warning: Sanitized betamax and betamin are more than 180 degrees apart.\n\tThat's probably not what you intended. betamax: " << commargs[BETAMAX]*180.0/M_PI << ", betamin: " << commargs[BETAMIN]*180.0/M_PI << "\n\tAre you trying to use put a beta range including zero? Edit the source code for that..." << std::endl;
        }
        //Now the input should be sanitized.
        //determine number of ranges for px, qx
        //as sin(alpha) can be negative -> fabs
        unsigned int maxn,maxm;
        maxn=fabs(commargs[B1MAX]/(commargs[A1]*sin(commargs[ALPHA])));
        maxm=fabs(commargs[B2MAX]/(commargs[A1]*sin(commargs[ALPHA])));
        //solutions run from -maxn to maxn, and from -maxm to maxm.
        //Due to the ambiguity of asin, two solutions exist for each value of n and m
        //This means, that (2*maxn+1)*2 solutions exist, the same for m.
        anglerange pxranges[4*maxn+2];
        anglerange qxranges[4*maxm+2];

        //up to now it was more or less dull C code. Now comes the first real difference:
        //As asin changes sign together with its argument, one has to treat positive and negative n differently
        //first the special case n=0:
        pxranges[0].setlower(commargs[ALPHA]);
        pxranges[0].setupper(commargs[ALPHA]);
        pxranges[1].setlower(commargs[ALPHA] - M_PI);
        pxranges[1].setupper(commargs[ALPHA] - M_PI);
        //now the slightly more difficult case: n>0
        for(i=1;i<=maxn;i++){
            //while factoring out the asin doesn't improve performance much - it's only used twice, it improves readability, as the important thing in the formulas below
            //are the signs.
            //the fmax and fmin are there, because maxn was calculated using B1MAX. Witht B1MIN the argument of arcsine can very well be outside its defined range.
            double asina1b1min=asin(fmax(-1.0,fmin(1.0,i*commargs[A1]*sin(commargs[ALPHA])/commargs[B1MIN])));
            double asina1b1max=asin(i*commargs[A1]*sin(commargs[ALPHA])/commargs[B1MAX]);
            //we need to consider that sin(alpha) can be negative. In that case the sign of the asin will change as well.
            if(asina1b1min>=0)
            {
                pxranges[2*i].setlower(commargs[ALPHA] - asina1b1min);
                pxranges[2*i].setupper(commargs[ALPHA] - asina1b1max);
                pxranges[2*i+1].setlower(commargs[ALPHA] - M_PI + asina1b1max);
                pxranges[2*i+1].setupper(commargs[ALPHA] - M_PI + asina1b1min);
                //and the most difficult case: n<0
                //here the arcsin is negative. as i is positive, I'll just change the sign in front of the arcsin.
                pxranges[2*i+2*maxn].setlower(commargs[ALPHA] + asina1b1max);
                pxranges[2*i+2*maxn].setupper(commargs[ALPHA] + asina1b1min);
                pxranges[2*i+2*maxn+1].setlower(commargs[ALPHA] - M_PI - asina1b1min);
                pxranges[2*i+2*maxn+1].setupper(commargs[ALPHA] - M_PI - asina1b1max);
            }
            else
            {
                //just as above, but with upper and lower limits switched
                pxranges[2*i].setlower(commargs[ALPHA] - asina1b1max);
                pxranges[2*i].setupper(commargs[ALPHA] - asina1b1min);
                pxranges[2*i+1].setlower(commargs[ALPHA] - M_PI + asina1b1min);
                pxranges[2*i+1].setupper(commargs[ALPHA] - M_PI + asina1b1max);
                //n<0
                pxranges[2*i+2*maxn].setlower(commargs[ALPHA] + asina1b1min);
                pxranges[2*i+2*maxn].setupper(commargs[ALPHA] + asina1b1max);
                pxranges[2*i+2*maxn+1].setlower(commargs[ALPHA] - M_PI - asina1b1max);
                pxranges[2*i+2*maxn+1].setupper(commargs[ALPHA] - M_PI - asina1b1min);
            }
        }
        //same nonsense for qxranges
        //first the easy part: m=0;
        qxranges[0].setlower(commargs[ALPHA] - commargs[BETAMAX]);
        qxranges[0].setupper(commargs[ALPHA] - commargs[BETAMIN]);
        qxranges[1].setlower(commargs[ALPHA] - commargs[BETAMAX] - M_PI);
        qxranges[1].setupper(commargs[ALPHA] - commargs[BETAMIN] - M_PI);
        //now the slightly more difficult case: m>0;
        for(i=1;i<=maxm;i++){
            //also here, the asin values are factored out for improved readability.
            double asina1b2min=asin(fmax(fmin(i*commargs[A1]*sin(commargs[ALPHA])/commargs[B2MIN],1.0),-1.0));
            double asina1b2max=asin(i*commargs[A1]*sin(commargs[ALPHA])/commargs[B2MAX]);
            //same here: keep in mind that sin(alpha) can be negative:
            if(asina1b2min>=0)
            {
                qxranges[2*i].setlower(commargs[ALPHA] - commargs[BETAMAX] - asina1b2min);
                qxranges[2*i].setupper(commargs[ALPHA] - commargs[BETAMIN] - asina1b2max);
                qxranges[2*i+1].setlower(commargs[ALPHA] - commargs[BETAMAX] - M_PI + asina1b2max);
                qxranges[2*i+1].setupper(commargs[ALPHA] - commargs[BETAMIN] - M_PI + asina1b2min);
                //and last, but not leasst, the most difficult, m<0 - here the arcsin is negative;
                //as i is positive, I'll just change the sign in front of the arcsin.
                qxranges[2*i+2*maxm].setlower(commargs[ALPHA] - commargs[BETAMAX] + asina1b2max);
                qxranges[2*i+2*maxm].setupper(commargs[ALPHA] - commargs[BETAMIN] + asina1b2min);
                qxranges[2*i+1+2*maxm].setlower(commargs[ALPHA] - commargs[BETAMAX] - M_PI - asina1b2min);
                qxranges[2*i+1+2*maxm].setupper(commargs[ALPHA] - commargs[BETAMIN] - M_PI - asina1b2max);
            }
            else
            {
                qxranges[2*i].setlower(commargs[ALPHA] - commargs[BETAMAX] - asina1b2max);
                qxranges[2*i].setupper(commargs[ALPHA] - commargs[BETAMIN] - asina1b2min);
                qxranges[2*i+1].setlower(commargs[ALPHA] - commargs[BETAMAX] - M_PI + asina1b2min);
                qxranges[2*i+1].setupper(commargs[ALPHA] - commargs[BETAMIN] - M_PI + asina1b2max);
                //m<0
                qxranges[2*i+2*maxm].setlower(commargs[ALPHA] - commargs[BETAMAX] + asina1b2min);
                qxranges[2*i+2*maxm].setupper(commargs[ALPHA] - commargs[BETAMIN] + asina1b2max);
                qxranges[2*i+1+2*maxm].setlower(commargs[ALPHA] - commargs[BETAMAX] - M_PI - asina1b2max);
                qxranges[2*i+1+2*maxm].setupper(commargs[ALPHA] - commargs[BETAMIN] - M_PI - asina1b2min);
            }
        }
        //Now the real big change from the plain C version: Using the anglerange class to generate a list
        //of overlaps instead of just outputting all of them.

        std::vector<anglerange> xoverlaps; //keep this, as the overlap between this and yoverlaps gives coincident results
        for(i=0;i<4*maxn+2;i++)
        {
            for(j=0;j<4*maxm+2;j++)
            {
                anglerange tmp = pxranges[i].overlap(qxranges[j]);
                if(!(tmp.isempty()))
                {
                    xoverlaps.push_back(tmp);
                }
            }
        }

        //ok, same thing for qy, py:
        //I know that this is wasteful regarding RAM, but well, we're still talking about bytes, not about megabytes ;-)
        unsigned int maxo, maxp;
        maxo=fabs(commargs[B1MAX]/(commargs[A2]*sin(commargs[ALPHA])));
        maxp=fabs(commargs[B2MAX]/(commargs[A2]*sin(commargs[ALPHA])));
        anglerange qyranges[4*maxo+2];
        anglerange pyranges[4*maxp+2];

        //now let's start with qyranges. As previously we need to consider the "sign" of o,p, and sin(alpha)
        //first: o=0
        qyranges[0].setlower(0.0);
        qyranges[0].setupper(0.0);
        qyranges[1].setlower(M_PI);
        qyranges[1].setupper(M_PI);
        for(i=1;i<=maxo;i++)
        {
            //also here: factor out the asin for improved readability.
            double asina2b1min = asin(fmax(-1.0,fmin(1.0,i*commargs[A2]*sin(commargs[ALPHA])/commargs[B1MIN])));
            double asina2b1max = asin(i*commargs[A2]*sin(commargs[ALPHA])/commargs[B1MAX]);
            //is sin(alpha)>0?
            if(asina2b1max>=0)
            {
                //case: o>0
                qyranges[2*i].setlower(asina2b1max);
                qyranges[2*i].setupper(asina2b1min);
                qyranges[2*i+1].setlower(M_PI - asina2b1min);
                qyranges[2*i+1].setupper(M_PI - asina2b1max);
                //case: o<0
                qyranges[2*i+2*maxo].setlower( -asina2b1min);
                qyranges[2*i+2*maxo].setupper( -asina2b1max);
                qyranges[2*i+1+2*maxo].setlower(M_PI + asina2b1max);
                qyranges[2*i+1+2*maxo].setupper(M_PI + asina2b1min);
            }
            else
            {
                //case: o>0
                qyranges[2*i].setlower(asina2b1min);
                qyranges[2*i].setupper(asina2b1max);
                qyranges[2*i+1].setlower(M_PI - asina2b1max);
                qyranges[2*i+1].setupper(M_PI - asina2b1min);
                //case: o<0
                qyranges[2*i+2*maxo].setlower( -asina2b1max);
                qyranges[2*i+2*maxo].setupper( -asina2b1min);
                qyranges[2*i+1+2*maxo].setlower(M_PI + asina2b1min);
                qyranges[2*i+1+2*maxo].setupper(M_PI + asina2b1max);
            }
        }
        //that was too easy. Probably it's buggy as hell...
        //now to py
        pyranges[0].setlower( -commargs[BETAMAX]);
        pyranges[0].setupper( -commargs[BETAMIN]);
        pyranges[1].setlower(M_PI - commargs[BETAMAX]);
        pyranges[1].setupper(M_PI - commargs[BETAMIN]);
        for(i=1;i<=maxp;i++)
        {
            //and again: readability
            double asina2b2min = asin(fmax(-1.0,fmin(1.0,i*commargs[A2]*sin(commargs[ALPHA])/commargs[B2MIN])));
            double asina2b2max = asin(i*commargs[A2]*sin(commargs[ALPHA])/commargs[B2MAX]);
            if(asina2b2max>=0)
            {
                //case: p>0
                pyranges[2*i].setlower(asina2b2max - commargs[BETAMAX]);
                pyranges[2*i].setupper(asina2b2min - commargs[BETAMIN]);
                pyranges[2*i+1].setlower(M_PI - asina2b2min - commargs[BETAMAX]);
                pyranges[2*i+1].setupper(M_PI - asina2b2max - commargs[BETAMIN]);
                //case: p<0
                pyranges[2*i+2*maxp].setlower( -asina2b2min - commargs[BETAMAX]);
                pyranges[2*i+2*maxp].setupper( -asina2b2max - commargs[BETAMIN]);
                pyranges[2*i+1+2*maxp].setlower(M_PI + asina2b2max - commargs[BETAMAX]);
                pyranges[2*i+1+2*maxp].setupper(M_PI + asina2b2min - commargs[BETAMIN]);
            }
            else
            {
                //ok, here the asin is of opposite sign!
                //case p>0
                pyranges[2*i].setlower(asina2b2min - commargs[BETAMAX]);
                pyranges[2*i].setupper(asina2b2max - commargs[BETAMIN]);
                pyranges[2*i+1].setlower(M_PI - asina2b2max - commargs[BETAMAX]);
                pyranges[2*i+1].setupper(M_PI - asina2b2min - commargs[BETAMIN]);
                //case: p<0
                pyranges[2*i+2*maxp].setlower( -asina2b2max - commargs[BETAMAX]);
                pyranges[2*i+2*maxp].setupper( -asina2b2min - commargs[BETAMIN]);
                pyranges[2*i+1+2*maxp].setlower(M_PI + asina2b2min - commargs[BETAMAX]);
                pyranges[2*i+1+2*maxp].setupper(M_PI + asina2b2max - commargs[BETAMIN]);
            }
        }
        //99 bottles of bugs on the wall, 99 bottles of bugs. You get one down and fix it up, 99 bottles of bugs...
        //100 bottles of bugs on the wall, 100 bottles of bugs....

        std::vector<anglerange> yoverlaps;
        for(i=0;i<4*maxo+2;i++)
        {
            for(j=0;j<4*maxp+2;j++)
            {
                anglerange tmp = qyranges[i].overlap(pyranges[j]);
                if(!(tmp.isempty()))
                {
                    yoverlaps.push_back(tmp);
                }
            }
        }

        //to be a commensurate match, an angle has to be in both, x- and yoverlaps
        std::vector<anglerange> commensurate;
        for(i=0;i<xoverlaps.size();i++)
        {
            for(j=0;j<yoverlaps.size();j++)
            {
                anglerange tmp = xoverlaps[i].overlap(yoverlaps[j]);
                if(!(tmp.isempty()))
                {
                    commensurate.push_back(tmp);
                }
            }
        }

        cout << "Possible coincident lattice matches with px, qx:" << std::endl;
        for(i=0;i<xoverlaps.size();i++)
        {
            cout << (xoverlaps.at(i)).getlower().getval()*180/M_PI << " " << (xoverlaps.at(i)).getupper().getval()*180/M_PI << std::endl;
        }
        cout << "Possible coincident lattice matches with qy, py:" << std::endl;
        for(i=0;i<yoverlaps.size();i++)
        {
            cout << (yoverlaps.at(i)).getlower().getval()*180/M_PI << " " << (yoverlaps.at(i)).getupper().getval()*180/M_PI << std::endl;
        }
        cout << "Possible commensurate lattice matches:" << std::endl;
        for(i=0;i<commensurate.size();i++)
        {
            cout << (commensurate.at(i)).getlower().getval()*180/M_PI << " " << (commensurate.at(i)).getupper().getval()*180/M_PI << std::endl;
        }
    }
}

