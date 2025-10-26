#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "tracing.hpp"
#include <stdlib.h>
#include <algorithm>

void determine_3D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices)
{
    // TRACE();
    /*3D reaction*/
    double Dr1 {};
    if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-10) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
    }

    double Dr2;
    if (std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-10) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
    }

    double Rmax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius };

    double sep {};
    double R1 {};
    bool withinRmax = get_distance(biMolData.pro1Index, biMolData.pro2Index, biMolData.relIface1, biMolData.relIface2,
        rxnIndex, rateIndex, isStateChangeBackRxn, sep, R1, Rmax, complexList, forwardRxns[rxnIndex], moleculeList, membraneObject);
    //     if(biMolData.pro1Index == 64)
    // 	std::cout <<"IN DETERMINE 3D, TRACKING 64 and 65! "<<biMolData.pro1Index<<" Rmax "<<Rmax<<" is within Rmax? "<<withinRmax<<" partner: "<<biMolData.pro2Index<<" sep "<<R1<<std::endl;

    //     if(biMolData.pro2Index == 64)
    // 	std::cout <<"IN DETERMINE 3D, TRACKING 64 and 65! "<<biMolData.pro1Index<<" Rmax "<<Rmax<<" is within Rmax? "<<withinRmax<<" partner: "<<biMolData.pro2Index<<" sep "<<R1<<std::endl;
    //     if(biMolData.pro1Index == 65)
    // 	std::cout <<"IN DETERMINE 3D, TRACKING 64 and 65! "<<biMolData.pro1Index<<" Rmax "<<Rmax<<" is within Rmax? "<<withinRmax<<" partner: "<<biMolData.pro2Index<<" sep "<<R1<<std::endl;

    //     if(biMolData.pro2Index == 65)
    // 	std::cout <<"IN DETERMINE 3D, TRACKING 64 and 65! "<<biMolData.pro1Index<<" Rmax "<<Rmax<<" is within Rmax? "<<withinRmax<<" partner: "<<biMolData.pro2Index<<" sep "<<R1<<std::endl;

    if (withinRmax) {
        // in case the molecule dissociated
        moleculeList[biMolData.pro1Index].probvec.push_back(0);
        moleculeList[biMolData.pro2Index].probvec.push_back(0);
    }

    if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated
        && moleculeList[biMolData.pro2Index].trajStatus != TrajStatus::propagated) {
        /*This movestat check is if you allow just dissociated proteins to avoid
         * overlap*/
        if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
            /*Evaluate probability of reaction, with reweighting*/
            double ratio { forwardRxns[rxnIndex].bindRadius / R1 };
            if (sep < 0) {
                if (biMolData.com1Index != biMolData.com2Index) {
                    // std::cout << "*****************************************************\n"
                    //           << " WARNING AT ITERATION " << simItr << "\n";
                    // std::cout << "SEPARATION BETWEEN INTERFACE " << biMolData.relIface1 << " ON MOLECULE "
                    //           << biMolData.pro1Index << " AND INTERFACE " << biMolData.relIface2 << " ON MOLECULE "
                    //           << biMolData.pro2Index << " IS LESS THAN 0\n";
                    // std::cout << "separation: " << sep << " r1: " << R1 << " p1: " << biMolData.pro1Index
                    //           << " p2: " << biMolData.pro2Index << " it " << simItr << " i1: " << biMolData.absIface1
                    //           << " i2: " << biMolData.absIface2 << '\n';
                    // std::cout << "MOL1 COM: " << moleculeList[biMolData.pro1Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro1Index].freelist.size() << '\n';
                    // std::cout << "IFACE1: "
                    //           << moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord << '\n';
                    // std::cout << "MOL2 COM: " << moleculeList[biMolData.pro2Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro2Index].freelist.size() << '\n';
                    // std::cout << "IFACE2: "
                    //           << moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord << '\n';
                    // std::cout << "*****************************************************\n";
                }

                sep = 0;
                ratio = 1;
                R1 = forwardRxns[rxnIndex].bindRadius;
            }

            /*If one particle (lipid) is bound in the membrane, reaction prob
            is half due to flux across only top half of particle.
            */
            double kdiff { 4 * M_PI * biMolData.Dtot * forwardRxns[rxnIndex].bindRadius };
            double kact { forwardRxns[rxnIndex].rateList[rateIndex].rate };



            /////////////////////////////      Change binding kinetics ////////////////////////////////////////////////
            //std::cout << "present_rxn_inx" <<  rxnIndex <<std::endl;


            //forwardRxns
            
            {
                //can have other data types instead of int but must same datatype as item 
                //std::find(vec.begin(), vec.end(), item) != vec.end()

                if (std::find(params.rxn_AA.begin(), params.rxn_AA.end(), rxnIndex) != params.rxn_AA.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);  
                    kact = params.k_AA[index];
                    //std::cout << rxnIndex << std::endl; 
                    //std::cout << kact << std::endl; 
                    //std::cout << "AA" << std::endl;

                }else if (std::find(params.rxn_AB.begin(), params.rxn_AB.end(), rxnIndex) != params.rxn_AB.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);   //Simit 1ms    index 1s 
                    kact = params.k_AB[index];

                }else if (std::find(params.rxn_BB.begin(), params.rxn_BB.end(), rxnIndex) != params.rxn_BB.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);   //Simit 1ms    index 1s 
                    kact = params.k_BB[index];
                    //std::cout << rxnIndex << std::endl; 
                    //std::cout << kact << std::endl; 
                    //std::cout << "BB" << std::endl;

                }else if (std::find(params.rxn_AC.begin(), params.rxn_AC.end(), rxnIndex) != params.rxn_AC.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);   //Simit 1ms    index 1s 
                    kact = params.k_AC[index];
                
                }else if (std::find(params.rxn_BC.begin(), params.rxn_BC.end(), rxnIndex) != params.rxn_BC.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);   //Simit 1ms    index 1s 
                    kact = params.k_BC[index];
                    //std::cout << rxnIndex << std::endl; 
                    //std::cout << kact << std::endl; 
                    //std::cout << "BC" << std::endl;

                }else if (std::find(params.rxn_CC.begin(), params.rxn_CC.end(), rxnIndex) != params.rxn_CC.end() ) {
                    int index = floor(simItr*params.timeStep/1e6);   //Simit 1ms    index 1s 
                    kact = params.k_CC[index];
                    //std::cout << rxnIndex << std::endl; 
                    //std::cout << kact << std::endl; 
                    //std::cout << "CC" << std::endl;    
                }
            }     
                


                /*
                for (int i; i<params.t_series.size(); i++){
                if (abs(params.t_series[i]*10.0-simItr)<=108.2){
                    kact = kact * params.k_AA[i]/20.0;
                    std::cout << simItr << std::endl;
                    std::cout << kact << std::endl;
                    break;
                }
                */
            
            

            /*
            //// 3D Spatial
            // adjust ka according to the RNA distribution around the particle. Refer to origin point 

            Coord iface11 = moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord;

            Coord iface22 = moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord;

            double r1 = iface11.get_magnitude();

            double r2 = iface22.get_magnitude();

            double disttmp = (r1 + r2) / 2.0;

            double gauss_sigma = membraneObject.waterBox.x/8.0;

            kact = kact * exp(- disttmp*disttmp/2.0/pow(gauss_sigma,2.0) );
            */
            

            /*
            //// Pseudo-2D Spatial 
            double xMean = (moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.x + moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.x)/2.0;
            double yMean = (moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.y + moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.y)/2.0;
            double zMean = (moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.z + moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.z)/2.0 - (-0.5*membraneObject.waterBox.z);

            double disttmp = pow(xMean * xMean + yMean * yMean + zMean * zMean, 0.5);

            double gauss_sigma = membraneObject.waterBox.x/8.0;

            kact = kact * exp(- disttmp*disttmp/2.0/pow(gauss_sigma,2.0) );
            */
           
          


            
            
            // adjust ka with temporal programming 
            
            
            //double timepresent = simItr*params.timeStep;
            //kact = kact * ((sin(  simItr  *2*M_PI/(100*1000)  )  + 1 )*0.5);  //period=2pi/constant     200 (*100)


            
            
            //for (int i; i<params.t_series.size(); i++){
            //    if (abs( (params.t_series[i]*10.00) - simItr )<=108.2){
            //        kact = kact*params.k_AA[i];
            //        break;
            //    }
            //}

                        
            /*
            if ( simItr >100*3000 && simItr <100*4000){
                kact=kact* ( (-1)/(1+exp(- (simItr-100*3000)  ))  +1 );
            }

            if ( simItr >100*4000 ){
                kact=kact*(1)/(1+exp(- (simItr-100*4000)  ));
            }            
            

            */
            
            
  
            

            /////////////////////////////      Change binding kinetics ////////////////////////////////////////////////



            if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-10 || std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-10)
                kact *= 2.0;

            if (forwardRxns[rxnIndex].isSymmetric == true)
                kact *= 2.0; // for A(a)+A(a)->A(a!).A(a!) case

            double fact { 1.0 + kact / kdiff };
            double alpha { fact * sqrt(biMolData.Dtot) / forwardRxns[rxnIndex].bindRadius };

            /*Check whether this pair was previously in reaction zone and
            therefore if reweighting should apply.
            Only need to store reweighting values for one protein in the pair,
            but for each pair need to know the protein partner and interface partner
            */
            double currnorm { 1.0 };
            double p0_ratio { 1.0 };

            /*protein i is molTypeIndex wprot and j is wprot2*/
            int proA = biMolData.pro1Index;
            int ifaceA = biMolData.relIface1;
            int proB = biMolData.pro2Index;
            int ifaceB = biMolData.relIface2;
            if (biMolData.pro1Index > biMolData.pro2Index) {
                proA = biMolData.pro2Index;
                ifaceA = biMolData.relIface2;
                proB = biMolData.pro1Index;
                ifaceB = biMolData.relIface1;
            }

            /*Find out what the previous reweighting was for this pair, if they
            were stored from previous step (inside reaction zone).
            */
            double rxnProb {};
            if (biMolData.com1Index != biMolData.com2Index) {
                /*don't renormalize if their positions are fixed by being in the same
                 * complex!*/
                for (int s { 0 }; s < moleculeList[proA].prevlist.size(); ++s) {
                    if (moleculeList[proA].prevlist[s] == proB && moleculeList[proA].prevmyface[s] == ifaceA
                        && moleculeList[proA].prevpface[s] == ifaceB) {
                        p0_ratio
                            = pirr_pfree_ratio_psF(R1, moleculeList[proA].prevsep[s], params.timeStep, biMolData.Dtot,
                                forwardRxns[rxnIndex].bindRadius, alpha, moleculeList[proA].ps_prev[s], 1E-10);
                        currnorm = moleculeList[proA].prevnorm[s] * p0_ratio;
                        s = moleculeList[proA].prevlist.size();
                    }
                }
                // std::cout << "calculating probvec1 for pros " << pro1 << ' ' << pro2
                // <<
                // '\n';
                rxnProb = passocF(R1, params.timeStep, biMolData.Dtot, forwardRxns[rxnIndex].bindRadius, alpha, kact / (kact + kdiff));
            } else {
                /*Proteins are in the same complex!*/
                /*If these proteins are at contact (sep=0) or close to at contact
                (<bindrad) then they will associate
                */
                if (sep < forwardRxns[rxnIndex].bindRadius)
                    rxnProb = 1.0;
            }
            if (rxnProb > 1.000001) {
                // std::cerr << "WARNING: prob of reaction is: " << rxnProb << " > 1. Avoid this using a smaller time step." << std::endl;
                //exit(1);
            }
            if (rxnProb > 0.5) {
                // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
            }

            moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;
            moleculeList[biMolData.pro2Index].probvec.back() = moleculeList[biMolData.pro1Index].probvec.back();

            moleculeList[proA].currprevsep.push_back(R1);
            moleculeList[proA].currlist.push_back(proB);
            moleculeList[proA].currmyface.push_back(ifaceA);
            moleculeList[proA].currpface.push_back(ifaceB);
            moleculeList[proA].currprevnorm.push_back(currnorm);
            moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
        } // Within reaction zone
    } // did not just dissociate
}
