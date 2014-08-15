#ifndef SAC_MODEL_ELLIPSE2D_H
#define SAC_MODEL_ELLIPSE2D_H

#include <pcl/sample_consensus/sac_model.h>
#include <pcl/sample_consensus/model_types.h>

namespace pcl {
/** \brief SampleConsensusModelEllipse2D defines a model for 2D ellipse segmentation on the X-Y plane.
  *
  * The model coefficients are defined as:
  *   - \b center.x : the X coordinate of the circle's center
  *   - \b center.y : the Y coordinate of the circle's center
  *   - \b a        : major axis radius
  *   - \b b        : minor axis radius
  *   - \b theta    : tilt angle of the major axis to the X axis
  *
  * \author Felix Mauch
  * \ingroup sample_consensus
  */
template <typename PointT>
class SampleConsensusModelEllipse2D : public SampleConsensusModel<PointT>
{
public:
    using SampleConsensusModel<PointT>::input_;
    using SampleConsensusModel<PointT>::indices_;
    using SampleConsensusModel<PointT>::radius_min_;
    using SampleConsensusModel<PointT>::radius_max_;
    using SampleConsensusModel<PointT>::error_sqr_dists_;

    typedef typename SampleConsensusModel<PointT>::PointCloud PointCloud;
    typedef typename SampleConsensusModel<PointT>::PointCloudPtr PointCloudPtr;
    typedef typename SampleConsensusModel<PointT>::PointCloudConstPtr PointCloudConstPtr;

    typedef boost::shared_ptr<SampleConsensusModelEllipse2D> Ptr;

    /** \brief Constructor for base SampleConsensusModelEllipse2D.
      * \param[in] cloud the input point cloud dataset
      * \param[in] random if true set the random seed to the current time, else set to 12345 (default: false)
      */
    SampleConsensusModelEllipse2D (const PointCloudConstPtr &cloud, bool random = false)
        : SampleConsensusModel<PointT> (cloud, random), tmp_inliers_ ()
    {};

    /** \brief Constructor for base SampleConsensusModelEllipse2D.
      * \param[in] cloud the input point cloud dataset
      * \param[in] indices a vector of point indices to be used from \a cloud
      * \param[in] random if true set the random seed to the current time, else set to 12345 (default: false)
      */
    SampleConsensusModelEllipse2D (const PointCloudConstPtr &cloud,
                                   const std::vector<int> &indices,
                                   bool random = false)
        : SampleConsensusModel<PointT> (cloud, indices, random), tmp_inliers_ ()
    {};

    /** \brief Copy constructor.
      * \param[in] source the model to copy into this
      */
    SampleConsensusModelEllipse2D (const SampleConsensusModelEllipse2D &source) :
        SampleConsensusModel<PointT> (), tmp_inliers_ ()
    {
        *this = source;
    }

    /** \brief Empty destructor */
    virtual ~SampleConsensusModelEllipse2D () {}

    /** \brief Copy constructor.
      * \param[in] source the model to copy into this
      */
    inline SampleConsensusModelEllipse2D&
    operator = (const SampleConsensusModelEllipse2D &source)
    {
        SampleConsensusModel<PointT>::operator=(source);
        tmp_inliers_ = source.tmp_inliers_;
        return (*this);
    }

    /** \brief Check whether the given index samples can form a valid 2D ellipse model, compute the model coefficients
      * from these samples and store them in model_coefficients. The ellipse coefficients are: x, y, r1x, r1y, r2x, r2y.
      * \param[in] samples the point indices found as possible good candidates for creating a valid model
      * \param[out] model_coefficients the resultant model coefficients
      */
    bool
    computeModelCoefficients (const std::vector<int> &samples,
                                    Eigen::VectorXf &model_coefficients);

    /** \brief Recompute the model coefficients using the given inlier set
      * and return them to the user. Pure virtual.
      *
      * @note: these are the coefficients of the model after refinement
      * (e.g., after a least-squares optimization)
      *
      * \param[in] inliers the data inliers supporting the model
      * \param[in] model_coefficients the initial guess for the model coefficients
      * \param[out] optimized_coefficients the resultant recomputed coefficients after non-linear optimization
      */
    void
    optimizeModelCoefficients (const std::vector<int> &inliers,
                               const Eigen::VectorXf &model_coefficients,
                               Eigen::VectorXf &optimized_coefficients);

    /** \brief Compute all distances from the cloud data to a given model. Pure virtual.
      *
      * \param[in] model_coefficients the coefficients of a model that we need to compute distances to
      * \param[out] distances the resultant estimated distances
      */
    void
    getDistancesToModel (const Eigen::VectorXf &model_coefficients,
                         std::vector<double> &distances);

    /** \brief Select all the points which respect the given model
      * coefficients as inliers. Pure virtual.
      *
      * \param[in] model_coefficients the coefficients of a model that we need to compute distances to
      * \param[in] threshold a maximum admissible distance threshold for determining the inliers from
      * the outliers
      * \param[out] inliers the resultant model inliers
      */
    void
    selectWithinDistance (const Eigen::VectorXf &model_coefficients,
                          const double threshold,
                          std::vector<int> &inliers);

    /** \brief Count all the points which respect the given model
      * coefficients as inliers. Pure virtual.
      *
      * \param[in] model_coefficients the coefficients of a model that we need to
      * compute distances to
      * \param[in] threshold a maximum admissible distance threshold for
      * determining the inliers from the outliers
      * \return the resultant number of inliers
      */
    virtual int
    countWithinDistance (const Eigen::VectorXf &model_coefficients,
                         const double threshold);

    /** \brief Create a new point cloud with inliers projected onto the model. Pure virtual.
      * \param[in] inliers the data inliers that we want to project on the model
      * \param[in] model_coefficients the coefficients of a model
      * \param[out] projected_points the resultant projected points
      * \param[in] copy_data_fields set to true (default) if we want the \a
      * projected_points cloud to be an exact copy of the input dataset minus
      * the point projections on the plane model
      */
    void
    projectPoints (const std::vector<int> &inliers,
                   const Eigen::VectorXf &model_coefficients,
                   PointCloud &projected_points,
                   bool copy_data_fields = true);

    /** \brief Verify whether a subset of indices verifies a given set of
      * model coefficients. Pure virtual.
      *
      * \param[in] indices the data indices that need to be tested against the model
      * \param[in] model_coefficients the set of model coefficients
      * \param[in] threshold a maximum admissible distance threshold for
      * determining the inliers from the outliers
      */
    bool
    doSamplesVerifyModel (const std::set<int> &indices,
                          const Eigen::VectorXf &model_coefficients,
                          const double threshold);

    /** \brief Check whether a model is valid given the user constraints.
      * \param[in] model_coefficients the set of model coefficients
      */
    bool
    isModelValid (const Eigen::VectorXf &model_coefficients);

    /** \brief Check if a sample of indices results in a good sample of points
      * indices. Pure virtual.
      * \param[in] samples the resultant index samples
      */
    bool
    isSampleGood (const std::vector<int> &samples) const;

    /** \brief Return an unique id for each type of model employed. */
    inline SacModel
    getModelType () const { return (SACMODEL_ELLIPSE2D); }

private:
    void removeColumnFromMatrix(const Eigen::MatrixXf& original, Eigen::MatrixXf& target, unsigned int col_nr);

    /** \brief Temporary pointer to a list of given indices for optimizeModelCoefficients () */
    const std::vector<int> *tmp_inliers_;

#if defined BUILD_Maintainer && defined __GNUC__ && __GNUC__ == 4 && __GNUC_MINOR__ > 3
#pragma GCC diagnostic warning "-Weffc++"
#endif
};

} // namespace pcl

#ifdef PCL_NO_PRECOMPILE
#include <pcl/sample_consensus/impl/sac_model_ellipse2d.hpp>
#endif

#endif // SAC_MODEL_ELLIPSE2D_H
