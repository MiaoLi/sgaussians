/*
 * main.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: Kim Seungsu
 */
#include <stdio.h>
#include "GMRDynamics.h"
#include "CDDynamics.h"

GMRDynamics *GMR_ready;
Gaussians *GMR_dir;
Gaussians *GMR_mag;

void initGMR(void)
{
	GMR_ready = new GMRDynamics(3, 6, 1./240.,
			"/home/seungsu/data/models/throwing/GMM_ready_mu.txt",
			"/home/seungsu/data/models/throwing/GMM_ready_sigma.txt",
			"/home/seungsu/data/models/throwing/GMM_ready_prio.txt" );
	GMR_ready->initGMR(0, 2, 3, 5);

/*
	GMR_dir = new Gaussians(3, 3*2,
			"/home/seungsu/data/models/throwing/GMR_dir_mu.txt",
			"/home/seungsu/data/models/throwing/GMR_dir_sigma.txt",
			"/home/seungsu/data/models/throwing/GMR_dir_prio.txt" );
	GMR_dir->InitFastGMR(0, 2, 3, 5);


	GMR_mag = new Gaussians(3, 4,
			"/home/seungsu/data/models/throwing/GMR_mag_mu.txt",
			"/home/seungsu/data/models/throwing/GMR_mag_sigma.txt",
			"/home/seungsu/data/models/throwing/GMR_mag_prio.txt" );
	GMR_mag->InitFastGMR(0,2,3,3);
	*/
}

void regressionGMR_ready(void)
{
	Vector lTarget(3);
	Vector lState(3), lNextState(3);

	lTarget(0) = 0.01;
	lTarget(1) = 0.01;
	lTarget(2) = 0.01;


	lState(0)= -0.1445;
	lState(1)= -0.0303;
	lState(2)= 0.3139;

	GMR_ready->setTarget(lTarget);
	GMR_ready->setState(lState);

	for(int i=0; i<50; i++){
		lNextState = GMR_ready->getNextState();
		printf("%5.3f %5.3f %5.3f \n", lNextState(0), lNextState(1), lNextState(2) );
	}

}

void testCDDyn(void)
{
	CDDynamics *lCDDyn;
	CDDynamics *lCDDynL;
	Vector lPos(3), lVel(3), lAccel(3), lTarget(3);
	Vector lAccelLimits(3);
	Vector lVelLimits(3);

	lPos.Zero();
	lTarget.One();
	lTarget *= 0.99;

	double lWn = 60.0;
	lCDDyn = new CDDynamics(3, 1./500., lWn);
	lCDDyn->SetStateTarget(lPos, lTarget);
	lCDDynL = new CDDynamics(3, 1./500., lWn);
	lCDDynL->SetStateTarget(lPos, lTarget);

	lAccelLimits.One();
	lAccelLimits *= 20.0; //
	lCDDynL->SetAccelLimits(lAccelLimits);

	lVelLimits.One();
	lVelLimits *= DEG2RAD(100.0);
	lCDDynL->SetVelocityLimits(lVelLimits);


	FILE *fid;
	fid =fopen("./log.txt", "w+");
	for(int i=0; i<500; i++)
	{
		lCDDyn->Update();
		lCDDyn->GetState(lPos, lVel);
		lCDDyn->GetStateAccel(lAccel);
		fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf ", lPos(0), lPos(1), lPos(2), lVel(0), lVel(1), lVel(2), lAccel(0), lAccel(1), lAccel(2));

		lCDDynL->Update();
		lCDDynL->GetState(lPos, lVel);
		lCDDynL->GetStateAccel(lAccel);
		fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", lPos(0), lPos(1), lPos(2), lVel(0), lVel(1), lVel(2), lAccel(0), lAccel(1), lAccel(2));

	}
	fclose(fid);

}

int main(int argc, char** argv)
{

	//initGMR();
	//regressionGMR_ready();

	testCDDyn();

	return 0;


}
