
#include <set>

#include <gmsh.h>

int main(int argc, char **argv)
{
  
    gmsh::initialize();
    gmsh::model::add("cilinder_model");
    double centerX = 0.0; // X-координата центра верха
    double centerY = 0.0; // Y-координата центра верха
    double centerZup = 1.0; // Z - верх
    double centerZdown = 0.0; // Z - низ центра
    double radius = 1.0;  // Радиус круга верха
    
    double lc = 0.1;

    // Создание точек для нижнего круга
    int centerunderTag = gmsh::model::geo::addPoint(centerX, centerY, centerZdown, lc);
    int startunderTag = gmsh::model::geo::addPoint(centerX + radius, centerY, centerZdown, lc);
    int endunderTag = gmsh::model::geo::addPoint(centerX - radius, centerY, centerZdown, lc);
    
    int circleunderTag1 = gmsh::model::geo::addCircleArc(startunderTag, centerunderTag, endunderTag); //1
    int circleunderTag2 = gmsh::model::geo::addCircleArc(endunderTag, centerunderTag, startunderTag); //2
    
    int curveLoopunderTag = gmsh::model::geo::addCurveLoop({circleunderTag1, circleunderTag2});
    
    //Создаем верхний круг
    
    int centerupTag = gmsh::model::geo::addPoint(centerX, centerY, centerZup, lc);
    int startupTag = gmsh::model::geo::addPoint(centerX + radius, centerY, centerZup, lc);
    int endupTag = gmsh::model::geo::addPoint(centerX - radius, centerY, centerZup, lc);
    
    int circleupTag1 = gmsh::model::geo::addCircleArc(startupTag, centerupTag, endupTag); //3
    int circleupTag2 = gmsh::model::geo::addCircleArc(endupTag, centerupTag, startupTag); //4
    
    int curveLoopupTag = gmsh::model::geo::addCurveLoop({circleupTag1, circleupTag2});
    
    // Создаем линии для поверхности посередине
    int firstmedlineTag = gmsh::model::geo::addLine(startunderTag, startupTag); //5
    
    int secondmedlineTag = gmsh::model::geo::addLine(endunderTag, endupTag); //6
    
    int curveLoopMedTag1 = gmsh::model::geo::addCurveLoop({circleupTag1, -secondmedlineTag, -circleunderTag1, firstmedlineTag});
    
    int curveLoopMedTag2 = gmsh::model::geo::addCurveLoop({circleupTag2, -firstmedlineTag, -circleunderTag2, secondmedlineTag});

    // Создание поверхностей
    
    
    // Создаем нижний круг
    int surfaceunderTag = gmsh::model::geo::addPlaneSurface({curveLoopunderTag});
    //Создаем верхний круг
    int surfaceupTag = gmsh::model::geo::addPlaneSurface({curveLoopupTag});
    
    // Создаем половину середины
    int surfacemediumTag1 = gmsh::model::geo::addSurfaceFilling({curveLoopMedTag1});
    // И вторую половину
    int surfacemediumTag2 = gmsh::model::geo::addSurfaceFilling({curveLoopMedTag2});
    
    //Остался Объем
    int SurfaceLoopTag1 = gmsh::model::geo::addSurfaceLoop({surfaceunderTag, surfacemediumTag1, surfacemediumTag2, surfaceupTag});
    
    int ValueTag = gmsh::model::geo::addVolume({SurfaceLoopTag1});
    // Синхронизация геометрии
    gmsh::model::geo::synchronize();

    // Генерация 2D сетки
    gmsh::model::mesh::generate(3);

  
    gmsh::write("circle_model.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

  
    gmsh::finalize();

    return 0;
}



