#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void create_complex_propagation_vectors(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    // the diffusion inside the sphere is the same as inside a box system.
    // on sphere surface, we set D.z = 0.

    if (targCom.D.z < 1E-14 && membraneObject.isSphere == true) { // on sphere surface
        Coord targTrans = create_complex_propagation_vectors_on_sphere(params, targCom);
        // now setup translational motion vector
        targCom.trajTrans.x = targTrans.x;
        targCom.trajTrans.y = targTrans.y;
        targCom.trajTrans.z = targTrans.z;

        // now setup the rotation
        targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
        targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
        targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();

    } else { // box system or inside the sphere (not on the sphere)
        // Create Gaussian distributed random translational motion
        targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
        targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
        const double PI = std::atan(1.0)*4;
        targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV() -37.4*1e-6* params.timeStep;  //50*1e-6 nm/us =  timestep (us)
        //-   params.timeStep*targCom.D.z*1e-18/(1.38e-23*318) *    targCom.memberList.size() * (4/3*PI*500e-9*500e-9*500e-9) * 55.0 *9.8 *1e9;  // D/kBT*U(X)  U(X)=mg'=  # * (4/3*pi^3)*55/1000*9.8    *55 is delta lo *unit conversion from m to nm

        // create Gaussian distributed random rotational motion
        targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
        targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
        targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();
    }

    //determine RS3Dinput
    double RS3Dinput { 0.0 };

    for (auto& molIndex : targCom.memberList) {
        for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
            if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
                RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                break;
            }
        }
    }

    //  reflect the boundary and also check_if_span, but not necessary for complex on sphere surface
    reflect_traj_complex_rad_rot(params, moleculeList, targCom, membraneObject, RS3Dinput);

    // // Create Gaussian distributed random translational motion
    // targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    // targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    // targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV();
    // // create Gaussian distributed random rotational motion
    // targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
    // targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
    // targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();
    // std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);
    // //reflect_traj_complex_rad_rot(params, moleculeList, targCom, M, membraneObject);

    // //determine RS3Dinput
    // double RS3Dinput { 0.0 };

    // for (auto& molIndex : targCom.memberList) {
    //     for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
    //         if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
    //             RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
    //             break;
    //         }
    //     }
    // }

    // reflect_traj_complex_rad_rot_new(params, moleculeList, targCom, M, membraneObject, RS3Dinput);
}
