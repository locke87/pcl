#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_ellipse2d.h>
#include <Eigen/Eigenvalues>

#include <pcl/io/pcd_io.h>



//using namespace pcl;

int main (int argc, char ** argv) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  // Create a test cloud
  // populate our PointCloud with points
  int length = 360;
  int length2 = 50;
  double a = 2.3;
  double b = 1.2;
  for (size_t i = 0; i < length; i=i++) {
//      for (size_t j = 0; j < length2; j++) {
          pcl::PointXYZ point(0,0,0);
          point.x = /*(double)j/length2 * */ a * std::cos(i);
          point.y = /*(double)j/length2 * */ b * std::sin(i);
          point.z = 0;
          cloud->push_back(point);
//      }
  }

//  float x_vec[] = {1, -1, -3, 3, 1};
//  float y_vec[] = {2, 2, 0, 0, -2};
//  for (int i=0; i<5; i++) {
//    pcl::PointXYZ point(0,0,0);
//    point.x = x_vec[i];
//    point.y = y_vec[i];
//    point.z = 0;
//    cloud->push_back(point);
//  }




  std::cout << "Created Pointcloud with " << cloud->size() << " points." << std::endl;
  pcl::io::savePCDFile("ellipse_cloud.pcd", *cloud);

  std::vector<int> inliers;

  pcl::SampleConsensusModelEllipse2D<pcl::PointXYZ>::Ptr model(new pcl::SampleConsensusModelEllipse2D<pcl::PointXYZ>(cloud));
  pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model);
  ransac.setDistanceThreshold(.1);
  std::cout << "Computing model..." << std::endl;
  ransac.computeModel(0);
  std::cout << "Looking for inliers..." << std::endl;
  ransac.getInliers(inliers);
  std::cout << "Number of found inliers: " << inliers.size() << std::endl;

  Eigen::VectorXf coefficients;
  ransac.getModelCoefficients(coefficients);

  double origin_error = sqrt(coefficients[0]*coefficients[0] + coefficients[1]*coefficients[1]);
  double major_axis_error = fabs(coefficients[2] - a);
  double minor_axis_error = fabs(coefficients[3] - b);
  double tilt_error = fabs(coefficients[4]);

  std::cout << "Origin error: " << origin_error << std::endl <<
               "Major axis error: " << major_axis_error << std::endl <<
               "Minor axis error: " << minor_axis_error << std::endl <<
               "Tilt error: " << tilt_error << std::endl;

  return 0;
}
