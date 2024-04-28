#include <llas.hpp>

int main(int argc, char** argv) {
#if defined(_WIN64)
  const auto& data = llas::read("C:\\Users\\araka\\Projects\\data\\Palac_Moszna.las", false);
  // const auto& data = llas::read("C:\\Users\\araka\\Projects\\data\\G_Sw_Anny.las", false);
#else
  const auto& data = llas::read("/home/shinaraka/Projects/llas/sample/Palac_Moszna.las");
#endif

  if (data) {
    const size_t& nPoints = data->getNumPoints();
    std::cout << "nPoints: " << nPoints << std::endl;

    const auto& coords = data->getPointCoords(false);
    std::cout << "coords.size(): " << coords.size() << std::endl;

    const auto& colors = data->getPointColors();
    std::cout << "colors.size(): " << colors.size() << std::endl;
  }

  return 0;
}
