#include "../../ThermoPropCubicEos.hpp"

int main() {
auto T = 350.; // Isoterma com temperatura constante

//auto databasePath = "/home/paloma/Downloads/ThermoPropCubicEos-main/ThermoPropCubicEos/database/test.yml";
//auto databasePath = "/home/paloma/Documentos/ThermoPropCubicEos/database/test.yml";

auto databasePath = "/home/paloma/Documentos/ThermoPropCubicEos/database/test.yml";

//"/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml"; 
auto components = "CO2 C1"; // Componentes
std::vector<double> Tc, Pc, omega;

//auto EoSModel = CubicEOSModel::PengRobinson; // PR escolhido
 //auto EoSModel = CubicEOSModel::VanDerWaals; // vdW escolhido
 auto EoSModel = CubicEOSModel::SoaveRedlichKwong;  //srk escolhido

std::vector<double> z{0.5,0.5}; // Fração molar para mistura
//std::vector<double> z{0.95,0.05}; //95% CO2
//std::vector<double> z{0.05,0.95}; 
auto ncomp=z.size();
read_database(Tc, Pc, omega, databasePath, components);
//srk 350 co2
std::vector<double> V = {0.0290533138, 0.0115695241, 0.0072096252, 0.0052293338, 0.0040981692, 0.0033664168, 0.0028542814, 0.0024758011, 0.0021846994, 0.0019538527, 0.0017663112, 0.0016109374, 0.0014801117, 0.0013684449, 0.0012720177, 0.0011879119, 0.0011139095, 0.0010482948, 0.0009897200, 0.0009371113, 0.0008896028, 0.0008464883, 0.0008071866, 0.0007712145, 0.0007381671, 0.0007077028, 0.0006795312, 0.0006534038, 0.0006291069, 0.0006064555, 0.0005852889, 0.0005654664, 0.0005468650, 0.0005293758, 0.0005129028, 0.0004973608, 0.0004826736, 0.0004687733, 0.0004555991, 0.0004430962, 0.0004312153, 0.0004199116, 0.0004091447, 0.0003988779, 0.0003890778, 0.0003797138, 0.0003707581, 0.0003621850, 0.0003539711, 0.0003460949, 0.0003385365, 0.0003312775, 0.0003243010, 0.0003175914, 0.0003111342, 0.0003049157, 0.0002989236, 0.0002931462, 0.0002875725, 0.0002821925, 0.0002769967, 0.0002719762, 0.0002671228, 0.0002624285, 0.0002578863, 0.0002534892, 0.0002492307, 0.0002451049, 0.0002411060, 0.0002372286, 0.0002334676, 0.0002298183, 0.0002262762, 0.0002228368, 0.0002194962, 0.0002162505, 0.0002130961, 0.0002100294, 0.0002070473, 0.0002041465, 0.0002013241, 0.0001985773, 0.0001959034, 0.0001932998, 0.0001907640, 0.0001882937, 0.0001858867, 0.0001835409, 0.0001812541, 0.0001790245, 0.0001768502, 0.0001747293, 0.0001726603, 0.0001706414, 0.0001686710, 0.0001667477, 0.0001648699, 0.0001630364, 0.0001612457, 0.0001594966};

/*
if (V.size() != 100) {
    std::cerr << "Erro: O vetor V deve conter exatamente 100 valores.\\n";
    return 1;
}
*/
calculateIsotermaComp(EoSModel, Tc, Pc, omega, T, V, z, ncomp); //utiliza regra de mistura
return 0;
}
/**/
