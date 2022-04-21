

//compiling on standard linux:
// ./clean.sh; g++ -o sim main.cpp -lm  -lgsl -lgslcblas

#define DEBUG
using namespace std;

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <sys/time.h>
#include "MersenneTwister.h"
#include <fstream>
#include <stdio.h>
#include <string>
#include "variables.cpp"
#include "write_data.cpp"



// function prototypes

MTRand mtrand1; //random number generator
int main()
{

//we simulate "rep" times the experiments labelling of "clones_num" number of cells.
for(int mm=0;mm<rep;mm++)
{
    
    // reinitializing arrays
    for(int i=0;i<4;i++)
    {
        count_pers[i]=0;
        for(int j=0;j<T;j++)
        {
            surv[i][j]=0;
            monoclonal[j]=0;
            av_size[j][0]=0;
            av_size[j][1]=0;
        }
    }
    
    
for(int m=0;m<clones_num;m++)
{
    for(int i=0;i<N1;i++)
    {
        for(int j=0;j<N2;j++)
        {
            type[i][j]=0; //on part de cellules sans type
            type_av[i][j]=0; //on part de cellules sans type
        }
    }

    // we simulate at the beginning a single clonal induction, picking a random cell in the crypt

int choose1=floor(4*mtrand1());
int choose2=floor(N2*mtrand1());
    
    // we keep track of the initial position of the labelling event
count_pers[choose1]+=1.0;
    
    // we actually label it ("1" in the crypt array)
type[choose1][choose2]=1;
type_av[choose1][choose2]=1;

    // initializing time
t=0;
while(t<T)	// one Monte-Carlo step in simulations
{
	for(int i=0;i<N;i++) // looping over all cells
	{


        //we pick a random cell to check if it divides
        ch1=floor(N1*mtrand1());
        ch2=floor(N2*mtrand1());
        double temp=mtrand1();
        if(temp<divrate)// then it does divide
        {
                double temp2=mtrand1(); // we randomly pick to know if the division is vertical or horizontal
                if(temp2>0.5) // then it's vertical, we move all cells in the column up accordingly
                {
                  
                    for(int j=ch1+1;j<N1;j++)
                    {
                        type[j][ch2]=type_av[j-1][ch2];
                    }
                
                    for(int j=0;j<N1;j++)
                    {
                        type_av[j][ch2]=type[j][ch2];
                    }
                }
                
                
                else // then the division was horizontal
                {
        
                    double temp3=mtrand1(); //we pick a random number to decide if it divides left or right
                    if(temp3>0.5) // division towards the right, all cells in the right column above are moved up by 1.
                    {
                        for(int j=ch1+1;j<N1;j++)
                        {
                            type[j][int((ch2+1)%N2)]=type_av[j-1][int((ch2+1)%N2)];
                        }
                        type[ch1][int((ch2+1)%N2)]=type_av[ch1][ch2];
                        
                        for(int i=0;i<N1;i++)
                        {
                            for(int j=0;j<N2;j++)
                            {
                                type_av[i][j]=type[i][j];
                            }
                        }
                    }
                    
                    if(temp3<0.5) //division towards the left, all cells in the left column above are moved up by 1.
                    {
                        int ch2_new=ch2-1;
                        if(ch2_new<0){ch2_new+=N2;}

                        for(int j=ch1+1;j<N1;j++)
                        {
                            type[j][ch2_new]=type_av[j-1][ch2_new];
                        }
                        type[ch1][ch2_new]=type_av[ch1][ch2];
                        
                        for(int i=0;i<N1;i++)
                        {
                            for(int j=0;j<N2;j++)
                            {
                                type_av[i][j]=type[i][j];
                            }
                        }
                        
                    }
                }
        }
        
        
        
        
        
        
        // we know check whether cells randomly intercalate
        temp=mtrand1();
        if(temp<inter) //intercalation event
        {
                double temp3=mtrand1();
                if((temp3>0.5 or ch1==0) and (ch1!=N1-1)) //then it's an intercalation with the cell above (except for the very top cell where it's not possible - if we have picked the bottom cell "ch==0" then we always have to intercalate with the above)
                {
                    ch1_new=ch1+1;
                    
                    type[ch1_new][ch2]=type_av[ch1][ch2];
                    type[ch1][ch2]=type_av[ch1_new][ch2];
                    
                    type_av[ch1][ch2]=type[ch1][ch2];
                    type_av[ch1_new][ch2]=type[ch1_new][ch2];
                }
                else //then it's an intercalation with the cell below
                {

                    ch1_new=ch1-1;
                    
                    type[ch1_new][ch2]=type_av[ch1][ch2];
                    type[ch1][ch2]=type_av[ch1_new][ch2];
                    
                    type_av[ch1][ch2]=type[ch1][ch2];
                    type_av[ch1_new][ch2]=type[ch1_new][ch2];
                    
                }   
            
        }
	}

    
    // counting the number of labelled cells in the crypts (first 4 rows)
    counts=0;
	for(int i=0;i<4;i++)
	{
        for(int j=0;j<N2;j++)
        {
            if(type[i][j]==1){counts+=1;}
        }
	}

    // if there are still a clonal footprint in the crypt, we add it to the persistence and the average clone size
    if(counts>0)
    {
        surv[choose1][int(t)]+=1.0;
        av_size[int(t)][0]+=counts;
        av_size[int(t)][1]+=1;
    }
    
    // if the crypt is monoclonal, we also record it
    if(counts==4*N2)
    {
        monoclonal[int(t)]+=1.0/clones_num;
    }
    
    
    t=t+1;

}
//end of the time loop for that clone

}
    // we write the data for that repetition
    write_data();
}
    

	return 0;
		
}
