#include <llas.hpp>

int main(int argc, char** argv) {
  // const auto& data = llas::read("C:\\Users\\araka\\Projects\\data\\Palac_Moszna.las");
  const auto& data = llas::read("/home/shinaraka/Projects/llas/sample/G_Sw_Anny.laz");

  const size_t& nPoints = data->getNumPoints();
  std::cout << "nPoints: " << nPoints << std::endl;

  const auto& coords = data->getPointCoords(false);
  std::cout << "coords.size(): " << coords.size() << std::endl;

  const auto& colors = data->getPointColors();
  std::cout << "colors.size(): " << colors.size() << std::endl;

  return 0;
}
