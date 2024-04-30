#include <llas.hpp>

int main(int argc, char** argv) {
#if defined(_WIN64)
  const std::vector<std::string> TEST_FILE_PATHS = {
      "C:\\Users\\araka\\Projects\\data\\Palac_Moszna.las",
      "C:\\Users\\araka\\Projects\\data\\G_Sw_Anny.las",
      "C:\\Users\\araka\\Projects\\data\\Trimble_StSulpice-Cloud-50mm.las"};
#else
  const std::vector<std::string> TEST_FILE_PATHS = {
      "/home/araka/Projects/llas/sample/Palac_Moszna.las",
      "/home/araka/Projects/llas/sample/G_Sw_Anny.las",
      "/home/araka/Projects/llas/sample/Trimble_StSulpice-Cloud-50mm.las"};
#endif

  for (const auto& filePath : TEST_FILE_PATHS) {
    const auto& data = llas::read(filePath);

    if (data) {
      const size_t& nPoints = data->getNumPoints();
      std::cout << "nPoints: " << nPoints << std::endl;

      const auto& coords = data->getPointCoords();
      std::cout << "coords.size(): " << coords.size() << std::endl;

      const auto& colors = data->getPointColors();
      std::cout << "colors.size(): " << colors.size() << std::endl;
    }
  }

  return 0;
}
