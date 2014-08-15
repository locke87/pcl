#ifndef PCL_SAMPLE_CONSENSUS_IMPL_SAC_MODEL_ELLIPSE_H_
#define PCL_SAMPLE_CONSENSUS_IMPL_SAC_MODEL_ELLIPSE_H_

#include <pcl/sample_consensus/sac_model_ellipse2d.h>
#include <pcl/sample_consensus/eigen.h>
#include <pcl/common/concatenate.h>
#include <Eigen/Eigenvalues>

//////////////////////////////////////////////////////////////////////////
template <typename PointT> bool
pcl::SampleConsensusModelEllipse2D<PointT>::computeModelCoefficients(const std::vector<int> &samples, Eigen::VectorXf &model_coefficients)
{
//    std::cout << "Computing ellipse model coefficients..." << std::endl;
    // Need 5 samples
    if (samples.size () != 5)
    {
        PCL_ERROR ("[pcl::SampleConsensusModelEllipse2D<PointT>::computeModelCoefficients] Invalid set of samples given (%lu)!\n", samples.size ());
        return (false);
    }

    model_coefficients.resize (5);

    Eigen::Vector2d p0 (input_->points[samples[0]].x, input_->points[samples[0]].y);
    Eigen::Vector2d p1 (input_->points[samples[1]].x, input_->points[samples[1]].y);
    Eigen::Vector2d p2 (input_->points[samples[2]].x, input_->points[samples[2]].y);
    Eigen::Vector2d p3 (input_->points[samples[3]].x, input_->points[samples[3]].y);
    Eigen::Vector2d p4 (input_->points[samples[4]].x, input_->points[samples[4]].y);


    // First calculate coefficients of a conic section through 5 points
    Eigen::MatrixXf m (5,6);
    m << p0[0]*p0[0], p0[0]*p0[1], p0[1]*p0[1], p0[0], p0[1], 1,
        p1[0]*p1[0], p1[0]*p1[1], p1[1]*p1[1], p1[0], p1[1], 1,
        p2[0]*p2[0], p2[0]*p2[1], p2[1]*p2[1], p2[0], p2[1], 1,
        p3[0]*p3[0], p3[0]*p3[1], p3[1]*p3[1], p3[0], p3[1], 1,
        p4[0]*p4[0], p4[0]*p4[1], p4[1]*p4[1], p4[0], p4[1], 1;

    Eigen::VectorXf conic_coefficients(6);
    for (int i=0; i<6; i++) {
      int sign =  -2*(i%2) + 1;
      Eigen::MatrixXf part = Eigen::MatrixXf::Zero(5,5);
      removeColumnFromMatrix(m, part, i);
      conic_coefficients[i] = sign * part.determinant();
    }


    // Check if coefficients can represent an ellipse
    if( (conic_coefficients[1]*conic_coefficients[1] - 4*conic_coefficients[0]*conic_coefficients[2]) < 0 ) {
      // Find ellipse parameters from conic coefficients
      Eigen::MatrixXf M0(3,3);
      M0 << conic_coefficients[5], 0.5*conic_coefficients[3], 0.5*conic_coefficients[4],
            0.5*conic_coefficients[3], conic_coefficients[0], 0.5*conic_coefficients[1],
            0.5*conic_coefficients[4], 0.5*conic_coefficients[1], conic_coefficients[2];
      Eigen::MatrixXf M(2,2);
      M << conic_coefficients[0], 0.5*conic_coefficients[1],
           0.5*conic_coefficients[1], conic_coefficients[2];
      double det_M0 = M0.determinant();
      double det_M = M.determinant();

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> es(M);
      Eigen::VectorXf lambda = es.eigenvalues();

      // origin_x
      model_coefficients[0] = (conic_coefficients[1]*conic_coefficients[4] - 2*conic_coefficients[2]*conic_coefficients[3]) /
                              (4*conic_coefficients[0]*conic_coefficients[2]-conic_coefficients[1]*conic_coefficients[1]);
      // origin_y
      model_coefficients[1] = (conic_coefficients[1]*conic_coefficients[3] - 2*conic_coefficients[0]*conic_coefficients[4]) /
                              (4*conic_coefficients[0]*conic_coefficients[2]-conic_coefficients[1]*conic_coefficients[1]);
      // major axis length
      model_coefficients[2] = std::sqrt(-1*det_M0 / (det_M * lambda[1]));
      // minor axis length
      model_coefficients[3] = std::sqrt(-1*det_M0 / (det_M * lambda[0]));
      // tilt angle
      model_coefficients[4] = ( M_PI/2 - std::atan((conic_coefficients[0]-conic_coefficients[2])/conic_coefficients[1])) / 2;

      return (true);
    } else {
      return (false);
    }
}

template <typename PointT> void
pcl::SampleConsensusModelEllipse2D<PointT>::optimizeModelCoefficients(const std::vector<int>& inliers, const Eigen::VectorXf& model_coefficients, Eigen::VectorXf& optimized_coefficients)
{

}

template <typename PointT> void
pcl::SampleConsensusModelEllipse2D<PointT>::getDistancesToModel(const Eigen::VectorXf& model_coefficients, std::vector<double>& distances)
{
  // Check if the model is valid given the user constraints
  if (!isModelValid (model_coefficients))
  {
    distances.clear ();
    return;
  }
  distances.resize (indices_->size ());

  double ell_cos = std::cos(model_coefficients[4]);
  double ell_sin = std::sin(model_coefficients[4]);
  double a_squared = model_coefficients[2] * model_coefficients[2];
  double b_squared = model_coefficients[3] * model_coefficients[3];

  // Iterate through the 3d points and calculate the distances from them to the ellipse
  for (size_t i = 0; i < indices_->size (); ++i)
  {
    // Transform query point into ellipse coorinate system first.
    // Then the distance is the same as if the point was in the first quadrant in the ellipse coordinate system
    double dx = input_->points[(*indices_)[i]].x - model_coefficients[0];
    double dy = input_->points[(*indices_)[i]].y - model_coefficients[1];
    double u = fabsf(dx * ell_cos - dy * ell_sin);
    double v = fabsf(dx * ell_sin - dy * ell_cos);

    // Find angle of ellipse point which is next to query point.
    // This loop converges towards the correct angle
    double phi = -1;
    double phi_next = 0;
    while (phi_next - phi > 0.01) {
      phi = phi_next;
      phi_next = std::atan( ((a_squared - b_squared) * std::sin(phi)  + v*model_coefficients[3]) / (u * model_coefficients[2]) );
    }

    // Claculate distance
    distances[i] = fabsf( sqrtf( (u - model_coefficients[2] * std::cos(phi)) * (u - model_coefficients[2] * std::cos(phi)) +
                                 (v - model_coefficients[3] * std::sin(phi)) * (v - model_coefficients[3] * std::sin(phi)) ));
  }

}

template <typename PointT> void
pcl::SampleConsensusModelEllipse2D<PointT>::selectWithinDistance(const Eigen::VectorXf& model_coefficients, const double threshold, std::vector<int>& inliers)
{
  // Check if the model is valid given the user constraints
  if (!isModelValid (model_coefficients))
  {
    inliers.clear ();
    return;
  }
  int nr_p = 0;
  inliers.resize (indices_->size ());
  error_sqr_dists_.resize (indices_->size ());

  std::vector<double> distances;
  getDistancesToModel(model_coefficients, distances);

  for (size_t i = 0; i < distances.size (); ++i)
  {
    if (distances[i] < threshold)
    {
      // Returns the indices of the points whose distances are smaller than the threshold
      inliers[nr_p] = (*indices_)[i];
      error_sqr_dists_[nr_p] = static_cast<double> (distances[i]);
      ++nr_p;
    }
  }
  inliers.resize (nr_p);
  error_sqr_dists_.resize (nr_p);
}

template <typename PointT> int
pcl::SampleConsensusModelEllipse2D<PointT>::countWithinDistance(const Eigen::VectorXf& model_coefficients, const double threshold)
{
  // Check if the model is valid given the user constraints
  if (!isModelValid (model_coefficients))
  {
    return (0);
  }
  int nr_p = 0;
  error_sqr_dists_.resize (indices_->size ());

  std::vector<double> distances;
  getDistancesToModel(model_coefficients, distances);

  for (size_t i = 0; i < distances.size (); ++i)
  {
    if (distances[i] < threshold)
    {
      ++nr_p;
    }
  }
  return (nr_p);
}

template <typename PointT> void
pcl::SampleConsensusModelEllipse2D<PointT>::projectPoints(const std::vector<int>& inliers, const Eigen::VectorXf& model_coefficients, pcl::SampleConsensusModelEllipse2D<PointT>::PointCloud& projected_points, bool copy_data_fields)
{
  // implement me!
}

template <typename PointT> bool
pcl::SampleConsensusModelEllipse2D<PointT>::doSamplesVerifyModel(const std::set<int>& indices, const Eigen::VectorXf& model_coefficients, const double threshold)
{
  // implement me!
}

template <typename PointT> bool
pcl::SampleConsensusModelEllipse2D<PointT>::isModelValid(const Eigen::VectorXf& model_coefficients)
{
  // Needs a valid model coefficients
    if (model_coefficients.size () != 5)
    {
      PCL_ERROR ("[pcl::SampleConsensusModelCircle2D::isModelValid] Invalid number of model coefficients given (%lu)!\n", model_coefficients.size ());
      return (false);
    }

    if (radius_min_ != -std::numeric_limits<double>::max() && model_coefficients[3] < radius_min_)
      return (false);
    if (radius_max_ != std::numeric_limits<double>::max() && model_coefficients[2] > radius_max_)
      return (false);

    return (true);
}

template <typename PointT> bool
pcl::SampleConsensusModelEllipse2D<PointT>::isSampleGood(const std::vector<int>& samples) const
{
  // implement me!
 return true;
}

template <typename PointT> void
pcl::SampleConsensusModelEllipse2D<PointT>::removeColumnFromMatrix(const Eigen::MatrixXf& original, Eigen::MatrixXf& target, unsigned int col_nr)
{
  if (original.rows() != target.rows() || original.cols() != target.cols()+1)
  {
    PCL_ERROR ("[pcl::SampleConsensusModelEllipse2D<PointT>::removeColumnFromMatrix] Matrix dimensions do not match. Target needs to have same row number and one column less than original.\n");
  }
  unsigned int num_rows = original.rows();
  unsigned int num_cols = original.cols();

  if( col_nr < num_cols ) {
      target.block(0,0,num_rows,col_nr) << original.block(0, 0, num_rows, col_nr);
      target.block(0,col_nr,num_rows,num_cols-col_nr-1) << original.block(0, col_nr+1, num_rows, num_cols - col_nr -1);
  }
}

#define PCL_INSTANTIATE_SampleConsensusModelEllipse2D(T) template class PCL_EXPORTS pcl::SampleConsensusModelEllipse2D<T>;

#endif    // PCL_SAMPLE_CONSENSUS_IMPL_SAC_MODEL_ELLIPSE_H_
