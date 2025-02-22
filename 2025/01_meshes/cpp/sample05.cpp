

#include <set>

#include <gmsh.h>

int main(int argc, char **argv)
{
  
    gmsh::initialize();
    gmsh::model::add("circle_model");
    double centerX = 0.0; // X-координата центра
    double centerY = 0.0; // Y-координата центра
    double radius = 1.0;  // Радиус круга

    // Создание точек для круга
    double lc = 1e-2;
    int centerTag = gmsh::model::geo::addPoint(centerX, centerY, 0.0, lc); // Центр круга
    int startTag = gmsh::model::geo::addPoint(centerX + radius, centerY, 0.0, lc); // Начальная точка окружности
    int endTag = gmsh::model::geo::addPoint(centerX - radius, centerY, 0.0, lc);
    // Создание дуги окружности
    int circleTag1 = gmsh::model::geo::addCircleArc(startTag, centerTag, endTag);
    int circleTag2 = gmsh::model::geo::addCircleArc(endTag, centerTag, startTag);
    // Создание кривой (окружности)
    int curveLoopTag = gmsh::model::geo::addCurveLoop({circleTag1, circleTag2});

    // Создание поверхности (круга)
    int surfaceTag = gmsh::model::geo::addPlaneSurface({curveLoopTag});

    // Синхронизация геометрии
    gmsh::model::geo::synchronize();

    // Генерация 2D сетки
    gmsh::model::mesh::generate(2);

  
    gmsh::write("circle_model.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

  
    gmsh::finalize();

    return 0;
}


