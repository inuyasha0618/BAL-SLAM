#include <iostream>
#include <fstream>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <cstdio>

using namespace std;

class BALProblem {
public:
    ~BALProblem() {
        delete[] point_index_;
        delete[] camera_index_;
        delete[] observations_;
        delete[] params_;
    }

    bool LoadFile(const char* filename) {
        FILE* fptr = fopen(filename, "r");

        if (fptr == NULL) {
            return false;
        }

        FscanfOrDie(fptr, "%d", &num_cameras_);
        FscanfOrDie(fptr, "%d", &num_points_);
        FscanfOrDie(fptr, "%d", &num_observations_);

        cout << num_cameras_ << " " << num_points_ << " " <<num_observations_ << endl;
        camera_index_ = new int[num_observations_];
        point_index_ = new int[num_observations_];
        observations_ = new double[2 * num_observations_];

        for (int i = 0; i < num_observations_; i++) {
            FscanfOrDie(fptr, "%d", camera_index_ + i);
            FscanfOrDie(fptr, "%d", point_index_ + i);
            FscanfOrDie(fptr, "%lf", observations_ + 2 * i);
            FscanfOrDie(fptr, "%lf", observations_ + 2 * i + 1);
        }

        num_params_ = 9 * num_cameras_ + 3 * num_points_;
        params_ =  new double[num_params_];

        for (int i = 0; i < num_params_; i++) {
            FscanfOrDie(fptr, "%lf", params_ + i);
        }

        return true;
    }

    void WriteToPLYFile(char* filename)const;

    int num_observations() {
        return num_observations_;
    }

    int mapObservationToCameraIndex(int observationIndex) {
        return camera_index_[observationIndex];
    }

    int mapObservationToPointIndex(int observationIndex) {
        return point_index_[observationIndex];
    }

    double getObservationXByIndex(int observationIndex) {
        double* pt = observations_ + 2 *observationIndex;
        return *pt;
    }

    double getObservationYByIndex(int observationIndex) {
        double* pt = observations_ + 2 *observationIndex + 1;
        return *pt;
    }

    double* cameraParamStart() const{
        return params_;
    }

    double* pointParamStart() const{
        return params_ + 9 * num_cameras_;
    }

    double* findCameraParamByObservationIndex(int observationIndex) {
        return params_ + 9 * mapObservationToCameraIndex(observationIndex);
    }

    double* findPointParamByObservationIndex(int observationIndex) {
        return pointParamStart() + 3 * mapObservationToPointIndex(observationIndex);
    }

    double* findCameraParamByCameraIndex(int cameraIndex) const{
        return params_ + 9 * cameraIndex;
    }

    double* findPointParamByPointIndex(int pointIndex) const{
        return pointParamStart() + 3 * pointIndex;
    }

    int num_cameras() const{
        return num_cameras_;
    }

    int num_points() const{
        return num_points_;
    }

private:
    int num_cameras_;
    int num_points_;
    int num_observations_;
    int num_params_;

    int* point_index_;
    int* camera_index_;
    double* observations_;
    double* params_;

    template <typename T>
    void FscanfOrDie(FILE* fptr, const char* format, T* value) {
        int num_scanned = fscanf(fptr, format, value);
        if (num_scanned != 1) {
            cout << "Error occured when scan file!" << endl;
        }
    }

};

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::WriteToPLYFile(char* filename)const{
  std::ofstream of(filename);

  of<< "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << num_cameras_ + num_points_
    << '\n' << "property float x"
    << '\n' << "property float y"
    << '\n' << "property float z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.

    double center[3];
    double invRot[3];
    for(int i = 0; i < num_cameras(); ++i){
      const double* camera = findCameraParamByCameraIndex(i);
      invRot[0] = -camera[0];
      invRot[1] = -camera[1];
      invRot[2] = -camera[2];

      ceres::AngleAxisRotatePoint(invRot, camera + 3, center);
      center[0] = -center[0];
      center[1] = -center[1];
      center[2] = -center[2];

      of << center[0] << ' ' << center[1] << ' ' << center[2]
         << "0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.

    for(int i = 0; i < num_points(); ++i){
      const double* point = findPointParamByPointIndex(i);
      for(int j = 0; j < 3; ++j){
        of << point[j] << ' ';
      }
      of << "255 255 255\n";
    }
    of.close();
}

struct ProjectionError {
    ProjectionError(double x, double y): observed_x(x), observed_y(y) {}

    template <typename T>
    bool operator()(const T* const camera, const T* const point, T* error) const {
        T p[3];
        ceres::AngleAxisRotatePoint(camera, point, p);

        p[0] += camera[3];
        p[1] += camera[4];
        p[2] += camera[5];

        T xp = - p[0] / p[2];
        T yp = - p[1] / p[2];

        T r2 = xp * xp + yp * yp;

        const T& f = camera[6];
        const T& k1 = camera[7];
        const T& k2 = camera[8];

        T distortion = T(1.0) + r2 * (k1 + k2 * r2);

        T computed_x = f * distortion * xp;
        T computed_y = f * distortion * yp;

        error[0] = T(observed_x) - computed_x;
        error[1] = T(observed_y) - computed_y;

        return true;
    }

    double observed_x;
    double observed_y;

    static ceres::CostFunction* CreateCost(const double observed_x, const double observed_y) {
        return new ceres::AutoDiffCostFunction<ProjectionError, 2, 9, 3>(new ProjectionError(observed_x, observed_y));
    }
};

int main(int argc, char** argv)
{
    ceres::Problem problem;
    BALProblem bal;

    cout << argv[1] << endl;

    if (!bal.LoadFile(argv[1])) {
        std::cerr << "ERROR: Unable to open file " << argv[1] << endl;
    }

    bal.WriteToPLYFile(argv[2]);

    for (int i = 0;i < bal.num_observations(); i++) {
        ceres::CostFunction* cost_function = ProjectionError::CreateCost(bal.getObservationXByIndex(i), bal.getObservationYByIndex(i));
        problem.AddResidualBlock(cost_function, NULL, bal.findCameraParamByObservationIndex(i), bal.findPointParamByObservationIndex(i));
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    bal.WriteToPLYFile(argv[3]);

    cout << summary.FullReport() << endl;

    return 0;
}

