# Light LAS (LLAS)
A simple, header-only, las-format file reader.

## Features
- Header-only
- Only used STL libraries
- Compatible with v1.4 LAS format (PointDataRecordFormat: 0 to 4)

## Usage
You only have to include `include/llas.hpp` file.

## Example
```cpp
#include "llas.hpp"

int main(int argc, char** argv) {
    // Read .las file
    const auto lasData = llas::read("sample.las");

    if (lasData != nullptr) {
        // You can access to Public Header.
        const auto systemIdentifier = lasData->publicHeader.systemIdentifier;
        
        // You can easily get point coordinates as linear vector.
        const auto pointCoords = data->getPointCoords();

        // You can easily get point colors as linear vector.
        const auto pointColors = data->getPointColors();
    }
}
```