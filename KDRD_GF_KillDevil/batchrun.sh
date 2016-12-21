#!/bin/bash

#// linear
#if(method_flag==0)
#{
#if(mess_flg==1) cout << "\n Initializing linear in x with " << npts
#<< " points\n" << endl;
#dx=(xm-x0)/(npts-1);
#// filling array
#for(xarr[0]=x0, i=1; i<npts; i++) xarr[i]=xarr[i-1]+dx;
#}
#
#// log
#if(method_flag==1)
#{
#if(mess_flg==1) cout << "\n Initializing logarithmic in x with " << npts
#<< " points\n" << endl;
#dx=pow(xm/x0, 1.0/(double)(npts-1));
#// filling array
#for(xarr[0]=x0, i=1; i<npts; i++) xarr[i]=xarr[i-1]*dx;
#}

nkpts=100
dx=$(echo )
for i in $(seq 1 0.2 5)
do
./KDRD $i
done
