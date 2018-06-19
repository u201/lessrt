# coding: utf-8
from xml.dom import minidom
import codecs
from Loger import log
import numpy as np

# Generate Atmosphere xml file


def readParametersFromDartMaketFile(dart_maket_file):
    layer_depths = []  # layer thickness for each layer
    g_params = []  # g parameters of aerosol for each band
    ext_coeff_mol = []  # for each layer and for each band, [[band1, band2,...],[band1, band2,...],...]
    single_albedo_mol = []
    ext_coeff_aerosol = []
    single_albedo_aerosol = []


    f = open(dart_maket_file)
    all_lines = f.readlines()

    # read atmosphere geometry
    arr = all_lines[1].split("\t")
    top_of_MA = float(arr[1])
    top_of_HA = float(arr[2])
    num_of_layer_MA = int(arr[5])
    num_of_layer_HA = int(arr[6])

    for i in range(num_of_layer_MA):
        layer_depths.append(top_of_MA/float(num_of_layer_MA))
    for i in range(num_of_layer_HA):
        layer_depths.append((top_of_HA - top_of_MA) / float(num_of_layer_HA))

    is_g_started = False
    is_layer_started = False


    for line in all_lines:
        # read g parameter
        if line.startswith(" * Atmosphere general properties"):
            is_g_started = True
            is_layer_started = False
            continue
        if line.startswith(" * Optical properties"):
            is_layer_started = True
            is_g_started = False
            ext_coeff_mol = [[] for i in range(len(g_params))]
            single_albedo_mol = [[] for i in range(len(g_params))]
            ext_coeff_aerosol = [[] for i in range(len(g_params))]
            single_albedo_aerosol = [[] for i in range(len(g_params))]
            continue

        if is_g_started:
            arr = line.split("\t")
            g_params.append(float(arr[7]))

        if is_layer_started:
            arr = line.split("\t")
            band_index = int(arr[0])
            ext_coeff_mol[band_index].append(float(arr[2]))
            # ext_coeff_mol[band_index].append(0)
            single_albedo_mol[band_index].append(float(arr[3]))
            # ext_coeff_aerosol[band_index].append(float(arr[4]))
            ext_coeff_aerosol[band_index].append(0)
            single_albedo_aerosol[band_index].append(float(arr[5]))

    f.close()
    ext_coeff_mol = list(map(list, zip(*ext_coeff_mol)))
    single_albedo_mol = list(map(list, zip(*single_albedo_mol)))
    ext_coeff_aerosol = list(map(list, zip(*ext_coeff_aerosol)))
    single_albedo_aerosol = list(map(list, zip(*single_albedo_aerosol)))
    return layer_depths, np.array(g_params), np.array(ext_coeff_mol), np.array(single_albedo_mol), np.array(ext_coeff_aerosol), np.array(single_albedo_aerosol)


def weighted_albedo(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol):
    rows, cols = ext_coeff_mol.shape
    total_albedo = np.zeros((rows, cols))
    for row in range(rows):
        for col in range(cols):
            ext_mol = ext_coeff_mol[row][col]
            ext_aerosol = ext_coeff_aerosol[row][col]
            albedo_mol = single_albedo_mol[row][col]
            albedo_aerosl = single_albedo_aerosol[row][col]
            if (ext_mol + ext_aerosol)  == 0:
                total_albedo[row][col] = 0
            else:
                total_albedo[row][col] = (ext_mol*albedo_mol+ext_aerosol*albedo_aerosl)/(ext_mol + ext_aerosol)
    return total_albedo

def weights_for_phase_function(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol):
    rows, cols = ext_coeff_mol.shape
    weights = np.zeros((rows, 2))  # for mol and aerosol
    for row in range(rows):
        ext_mol = ext_coeff_mol[row].mean()
        ext_aerosol = ext_coeff_aerosol[row].mean()
        albedo_mol = single_albedo_mol[row].mean()
        albedo_aerosl = single_albedo_aerosol[row].mean()
        total_weights = ext_mol*albedo_mol+ext_aerosol*albedo_aerosl
        if total_weights == 0:
            weights[row][0] = 0.5
            weights[row][1] = 0.5
        else:
            weights[row][0] = (ext_mol*albedo_mol)/(ext_mol*albedo_mol+ext_aerosol*albedo_aerosl)
            weights[row][1] = (ext_aerosol * albedo_aerosl) / (ext_mol * albedo_mol + ext_aerosol * albedo_aerosl)
    return weights


def generateAtsXmlNF(dart_maket_path, output_xml_path, width, height):
    log("INFO: Generating atmosphere.")
    layer_depths, g_params, ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol \
        = readParametersFromDartMaketFile(dart_maket_path)

    total_ext_coeff = ext_coeff_mol + ext_coeff_aerosol
    total_albedo = weighted_albedo(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol)
    weightsForPhase = weights_for_phase_function(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol,
                                                 single_albedo_aerosol)

    f = codecs.open(output_xml_path, "w", "utf-8-sig")
    doc = minidom.Document()
    root = doc.createElement("scene")
    doc.appendChild(root)
    root.setAttribute("version", "0.5.0")

    vertical_offset = -0.00001

    # shapeNode = doc.createElement("shape")
    # root.appendChild(shapeNode)
    # shapeNode.setAttribute("type","cube")
    # toWorldNode = doc.createElement("transform")
    # shapeNode.appendChild(toWorldNode)
    # toWorldNode.setAttribute("name", "toWorld")
    # scaleNode = doc.createElement("scale")
    # toWorldNode.appendChild(scaleNode)
    # scaleNode.setAttribute("x", str(width * 0.5))
    # scaleNode.setAttribute("z", str(height * 0.5))
    # scaleNode.setAttribute("y", str(sum(layer_depths) * 0.5))  # total atmosphere height
    # translateNode = doc.createElement("translate")
    # toWorldNode.appendChild(translateNode)
    # translateNode.setAttribute("y", str(sum(layer_depths) * 0.5 + vertical_offset))

    mediumNode = doc.createElement("medium")
    root.appendChild(mediumNode)
    mediumNode.setAttribute("id","ats_medium")
    mediumNode.setAttribute("name", "interior")
    mediumNode.setAttribute("type", "planeparallel")
    floatNode = doc.createElement("float")
    mediumNode.appendChild(floatNode)
    floatNode.setAttribute("name","startAltitude")
    floatNode.setAttribute("value",str(vertical_offset))
    strNode = doc.createElement("string")
    mediumNode.appendChild(strNode)
    strNode.setAttribute("name","layerThickness")
    strNode.setAttribute("value",",".join([str(t) for t in layer_depths]))
    layer_index = 1
    for i in range(len(layer_depths)):
        if any(total_ext_coeff[i]>0):
            specNode = doc.createElement("spectrum")
            mediumNode.appendChild(specNode)
            specNode.setAttribute("name","singmaT_layer"+str(layer_index))
            specNode.setAttribute("value",",".join([str(t) for t in total_ext_coeff[i]]))
            specNode = doc.createElement("spectrum")
            mediumNode.appendChild(specNode)
            specNode.setAttribute("name", "albedo_layer" + str(layer_index))
            specNode.setAttribute("value", ",".join([str(t) for t in total_albedo[i]]))
            strNode = doc.createElement("string")
            mediumNode.appendChild(strNode)
            strNode.setAttribute("name","phasefunc_weights_layer"+str(layer_index))
            strNode.setAttribute("value",",".join([str(t) for t in weightsForPhase[i]]))
            floatNode = doc.createElement("float")
            mediumNode.appendChild(floatNode)
            floatNode.setAttribute("name","phasefunc_g_value_layer"+str(layer_index))
            floatNode.setAttribute("value",str(g_params.mean()))
            layer_index += 1


    xm = doc.toprettyxml()
    # xm = xm.replace('<?xml version="1.0" ?>', '')
    f.write(xm)
    f.close()
    log("INFO: Atmosphere generated.")


def generateAtsXml(dart_maket_path, output_xml_path, width, height):
    log("INFO: Generating atmosphere.")
    layer_depths, g_params, ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol\
        = readParametersFromDartMaketFile(dart_maket_path)

    total_ext_coeff = ext_coeff_mol + ext_coeff_aerosol
    total_albedo = weighted_albedo(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol)
    weightsForPhase = weights_for_phase_function(ext_coeff_mol, single_albedo_mol, ext_coeff_aerosol, single_albedo_aerosol)


    f = codecs.open(output_xml_path, "w", "utf-8-sig")
    doc = minidom.Document()
    root = doc.createElement("scene")
    doc.appendChild(root)
    root.setAttribute("version", "0.5.0")

    layer_bottom_height = [0]
    vertical_offset = 2
    for i in range(1, len(layer_depths)):
        bottom_height = sum(layer_depths[0:i])
        layer_bottom_height.append(bottom_height)
    for i in range(len(layer_depths)):
        shapeNode = doc.createElement("shape")
        root.appendChild(shapeNode)
        shapeNode.setAttribute("type","cube")
        toWorldNode = doc.createElement("transform")
        shapeNode.appendChild(toWorldNode)
        toWorldNode.setAttribute("name","toWorld")
        scaleNode = doc.createElement("scale")
        toWorldNode.appendChild(scaleNode)
        scaleNode.setAttribute("x", str(width*0.5))
        scaleNode.setAttribute("z",str(height*0.5))
        scaleNode.setAttribute("y", str(layer_depths[i] * 0.5-0.1))
        translateNode  = doc.createElement("translate")
        toWorldNode.appendChild(translateNode)
        translateNode.setAttribute("y", str(layer_bottom_height[i]+layer_depths[i] * 0.5+vertical_offset))
        mediumNode = doc.createElement("medium")
        shapeNode.appendChild(mediumNode)
        mediumNode.setAttribute("name","interior")
        mediumNode.setAttribute("type","homogeneous")
        spectrumNode = doc.createElement("spectrum")
        mediumNode.appendChild(spectrumNode)
        spectrumNode.setAttribute("name","sigmaT")
        spectrumNode.setAttribute("value",",".join([str(t) for t in total_ext_coeff[i]]))
        spectrumNode = doc.createElement("spectrum")
        mediumNode.appendChild(spectrumNode)
        spectrumNode.setAttribute("name","albedo")
        spectrumNode.setAttribute("value",",".join([str(t) for t in total_albedo[i]]))

        phaseNode = doc.createElement("phase")
        mediumNode.appendChild(phaseNode)
        phaseNode.setAttribute("type","mixturephase")
        weightsNode = doc.createElement("string")
        phaseNode.appendChild(weightsNode)
        weightsNode.setAttribute("name","weights")
        weightsNode.setAttribute("value",",".join([str(t) for t in weightsForPhase[i]]))
        rayleighNode = doc.createElement("phase")
        phaseNode.appendChild(rayleighNode)
        rayleighNode.setAttribute("type","rayleigh")
        hgNode = doc.createElement("phase")
        phaseNode.appendChild(hgNode)
        hgNode.setAttribute("type","hg")
        floatNode = doc.createElement("float")
        hgNode.appendChild(floatNode)
        floatNode.setAttribute("name","g")
        floatNode.setAttribute("value",str(g_params.mean()))

    xm = doc.toprettyxml()
    # xm = xm.replace('<?xml version="1.0" ?>', '')
    f.write(xm)
    f.close()
    log("INFO: Atmosphere generated.")

if __name__ == "__main__":
    readParametersFromDartMaketFile(r'D:\DART\user_data\simulations\AtsForLess\output\atmosphereMaket.txt')