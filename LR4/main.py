from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from model_export import *
from reflection_models import *
import numpy as np
import json
from PIL import Image  # Для сохранения BMP
import sys

window_width = 1020
window_height = 820


def read_json(json_file):
    with open(json_file, 'r') as file:
        j = json.load(file)
    return j['name'], j


def read_DeapthMap(depthmap_file):
    with open(depthmap_file, 'rb') as file:
        height = int(np.fromfile(file, dtype=np.float64, count=1)[0])
        width = int(np.fromfile(file, dtype=np.float64, count=1)[0])
        depthmap_array = np.fromfile(
            file, dtype=np.float64).reshape((height, width))
    return depthmap_array

# вычисление нормали
def calculate_normal(p1, p2, p3):
    u = np.array(p2) - np.array(p1)
    v = np.array(p3) - np.array(p1)
    normal = np.cross(u, v)

    norm = np.linalg.norm(normal)
    if norm == 0:
        return normal

    return normal / norm

# выбор формата экспорта
def export(depthmap_array, format, name):
    output_path = f"{name}.{format}"
    if format == 'vrml':
        export_vrml(output_path, depthmap_array)
    elif format == 'stl':
        export_stl(output_path, depthmap_array)
    elif format == 'ply':
        export_ply(output_path, depthmap_array)

# экспорт в формат VRML
def export_vrml(file, depthmap_array):
    height, width = depthmap_array.shape
    vertices = []
    indices = {}
    faces = []

    # вершины
    for i in range(height):
        for j in range(width):
            depth = depthmap_array[i, j]
            if depth != 0:
                vertex = (j - width / 2, height / 2 - i, -depth)
                indices[(i, j)] = len(vertices)
                vertices.append(vertex)

    # грани
    for i in range(height - 1):
        for j in range(width - 1):
            if (depthmap_array[i, j] != 0 and
                    depthmap_array[i, j + 1] != 0 and
                    depthmap_array[i + 1, j] != 0 and
                    depthmap_array[i + 1, j + 1] != 0):
                v1 = indices[(i, j)]
                v2 = indices[(i, j + 1)]
                v3 = indices[(i + 1, j + 1)]
                v4 = indices[(i + 1, j)]
                faces.append((v1, v2, v3, v4))

    with open(file, 'w') as file:
        file.write("#VRML V2.0 utf8\n\n")
        file.write("Shape {\n")
        file.write("  appearance Appearance {\n")
        file.write("    material Material {\n")
        file.write("      diffuseColor 0.8 0.8 0.8\n")
        file.write("      specularColor 0.5 0.5 0.5\n")
        file.write("      shininess 0.5\n")
        file.write("    }\n")
        file.write("  }\n")
        file.write("  geometry IndexedFaceSet {\n")
        file.write("    solid FALSE\n")
        file.write("    coord Coordinate {\n")
        file.write("      point [\n")

        for vertex in vertices:
            file.write(f"        {vertex[0]:.6f} {vertex[1]:.6f} {vertex[2]:.6f},\n")

        file.write("      ]\n")
        file.write("    }\n")
        file.write("    coordIndex [\n")

        for face in faces:
            file.write(f"        {face[0]}, {face[1]}, {face[2]}, {face[3]}, -1,\n")

        file.write("    ]\n")
        file.write("  }\n")
        file.write("}\n")

        print(f"  VRML: {len(vertices)} вершин, {len(faces)} граней")
# экспорт в формат STL
def export_stl(file, depthmap_array):
    height, width = depthmap_array.shape
    with open(file, 'w') as file:
        file.write('solid depth_map\n')
        for i in range(height - 1):
            for j in range(width - 1):
                if (depthmap_array[i, j] != 0 and
                        depthmap_array[i, j + 1] != 0 and
                        depthmap_array[i + 1, j] != 0 and
                        depthmap_array[i + 1, j + 1] != 0):

                    vertices = [
                        [j - width / 2, height / 2 - i, -depthmap_array[i, j]],
                        [j + 1 - width / 2, height / 2 - i, -depthmap_array[i, j + 1]],
                        [j + 1 - width / 2, height / 2 - (i + 1), -depthmap_array[i + 1, j + 1]],
                        [j - width / 2, height / 2 - (i + 1), -depthmap_array[i + 1, j]]
                    ]

                    for k in range(2):
                        triangle = (vertices[0], vertices[k + 1], vertices[k + 2])
                        normal = calculate_normal(*triangle)
                        normal_string = ' '.join(map(str, normal))
                        file.write(f'facet normal {normal_string}\n')
                        file.write(' outer loop\n')
                        for vertex in triangle:
                            vertex_string = ' '.join(map(str, vertex))
                            file.write(f' vertex {vertex_string}\n')
                        file.write(' endloop\n')
                        file.write('endfacet\n')
        file.write('endsolid depth_map\n')

# экспорт в формат ply
def export_ply(file, depthmap_array):
    height, width = depthmap_array.shape
    vertices = []
    indices = {}
    faces = []

    for i in range(height):
        for j in range(width):
            depth = depthmap_array[i, j]
            if depth != 0:
                vertex = (j - width / 2, height / 2 - i, -depth)
                indices[(i, j)] = len(vertices)
                vertices.append(vertex)

    for i in range(height - 1):
        for j in range(width - 1):
            if (depthmap_array[i, j] != 0 and
                    depthmap_array[i, j + 1] != 0 and
                    depthmap_array[i + 1, j] != 0 and
                    depthmap_array[i + 1, j + 1] != 0):
                v1 = indices[(i, j)]
                v2 = indices[(i, j + 1)]
                v3 = indices[(i + 1, j + 1)]
                v4 = indices[(i + 1, j)]

                faces.append((v1, v2, v3, v4))

    with open(file, 'w') as file:
        file.write("ply\n")
        file.write("format ascii 1.0\n")
        file.write(f"element vertex {len(vertices)}\n")
        file.write("property float x\n")
        file.write("property float y\n")
        file.write("property float z\n")
        file.write(f"element face {len(faces)}\n")
        file.write("property list uchar int vertex_index\n")
        file.write("end_header\n")

        for vertex in vertices:
            file.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")

        for face in faces:
            file.write(f"4 {face[0]} {face[1]} {face[2]} {face[3]}\n")

# выбор модель отражения
def models(normal, light, camera, model):
    if model == 'lambert':
        return lambert(normal, light)
    elif model == 'phong':
        return phong(normal, light, camera)
    elif model == 'oren-nayar':
        return oren_nayar(normal, light, camera)
    else:
        return lambert(normal, light)

# модель Ламберта
def lambert(normal, light_source):
    light_vector = np.array(light_source) / np.linalg.norm(light_source)
    normal_vector = np.array(normal) / np.linalg.norm(normal)
    dot = np.dot(normal_vector, light_vector)
    intensity = max(dot, 0)
    intensity = 0.1 + 0.9 * intensity
    intensity = min(intensity, 1)
    return (intensity, intensity, intensity)

# модель Орена-Найара
def oren_nayar(normal, light_source, viewer):
    # Нормализация векторов
    light_vector = np.array(light_source) / np.linalg.norm(light_source)
    normal_vector = np.array(normal) / np.linalg.norm(normal)
    viewer_vector = np.array(viewer) / np.linalg.norm(viewer)

    # Параметры модели Орена-Найара
    roughness = 0.8  # Коэффициент шероховатости (0-1)
    albedo = 0.7  # Альбедо поверхности

    # Углы падения и отражения
    cos_theta_i = max(0.0, min(1.0, np.dot(normal_vector, light_vector)))
    cos_theta_r = max(0.0, min(1.0, np.dot(normal_vector, viewer_vector)))

    theta_i = math.acos(cos_theta_i)
    theta_r = math.acos(cos_theta_r)

    # Проекции на плоскость, перпендикулярную нормали
    light_projection = light_vector - normal_vector * cos_theta_i
    viewer_projection = viewer_vector - normal_vector * cos_theta_r

    # Косинус разности азимутальных углов
    if np.linalg.norm(light_projection) > 1e-10 and np.linalg.norm(viewer_projection) > 1e-10:
        light_projection = light_projection / np.linalg.norm(light_projection)
        viewer_projection = viewer_projection / np.linalg.norm(viewer_projection)
        cos_phi_diff = np.dot(light_projection, viewer_projection)
        cos_phi_diff = max(-1.0, min(1.0, cos_phi_diff))
    else:
        cos_phi_diff = 0.0

    # Параметры модели
    sigma = roughness
    sigma2 = sigma * sigma

    A = 1.0 - 0.5 * (sigma2 / (sigma2 + 0.33))
    B = 0.45 * (sigma2 / (sigma2 + 0.09))

    alpha = max(theta_i, theta_r)
    beta = min(theta_i, theta_r)

    # Основная формула Орена-Найара
    term1 = albedo * cos_theta_i
    term2 = A + B * max(0.0, cos_phi_diff) * math.sin(alpha) * math.tan(beta)

    intensity = term1 * term2

    # Гарантируем минимальную видимость
    intensity = 0.1 + 0.9 * max(0.0, intensity)
    intensity = min(intensity, 1.0)

    return (intensity, intensity, intensity)

# модель Фонга
def phong(normal, light_source, viewer):
    light_vector = np.array(light_source) / np.linalg.norm(light_source)
    normal_vector = np.array(normal) / np.linalg.norm(normal)
    viewer_vector = np.array(viewer) / np.linalg.norm(viewer)

    ambient = 0.1
    diffuse_coefficient = 0.7
    specular_coefficient = 0.2
    shininess = 32

    ambient_color = ambient

    dot_l_n = max(np.dot(light_vector, normal_vector), 0)
    diffuse_color = diffuse_coefficient * dot_l_n

    reflection_vector = 2 * normal_vector * dot_l_n - light_vector
    dot_r_v = max(np.dot(reflection_vector, viewer_vector), 0)
    specular_color = specular_coefficient * (dot_r_v ** shininess)

    intensity = ambient_color + diffuse_color + specular_color
    intensity = min(intensity, 1)

    return (intensity, intensity, intensity)

def render_scene(depthmap_array, light, camera, shading):
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    gluLookAt(camera[0], camera[1], camera[2],
              0, 0, 0,
              0, 1, 0)

    height, width = depthmap_array.shape
    max_depth = np.max(depthmap_array) if np.max(depthmap_array) != 0 else 1
    scale_factor = 200

    glBegin(GL_QUADS)
    for i in range(height - 1):
        for j in range(width - 1):
            if (depthmap_array[i, j] != 0 and
                    depthmap_array[i, j + 1] != 0 and
                    depthmap_array[i + 1, j] != 0 and
                    depthmap_array[i + 1, j + 1] != 0):

                vertices = [
                    ((j - width / 2) / width, (height / 2 - i) / height, -depthmap_array[i, j] / max_depth),
                    (((j + 1) - width / 2) / width, (height / 2 - i) / height, -depthmap_array[i, j + 1] / max_depth),
                    (((j + 1) - width / 2) / width, (height / 2 - (i + 1)) / height,
                     -depthmap_array[i + 1, j + 1] / max_depth),
                    ((j - width / 2) / width, (height / 2 - (i + 1)) / height, -depthmap_array[i + 1, j] / max_depth)
                ]

                normal = calculate_normal(vertices[0], vertices[1], vertices[2])
                color = models(normal, light, camera, shading)

                glColor3f(*color)
                for vertex in vertices:
                    glVertex3f(vertex[0] * scale_factor, vertex[1] * scale_factor, vertex[2] * scale_factor)
    glEnd()
    glutSwapBuffers()


def init_glut(display_func):
    glutInit()
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowSize(window_width, window_height)
    glutInitWindowPosition(100, 100)
    glutCreateWindow(b"3D")
    glutDisplayFunc(display_func)
    glEnable(GL_DEPTH_TEST)
    glClearColor(0.0, 0.0, 0.0, 0.0)

    glMatrixMode(GL_PROJECTION)
    gluPerspective(45, (window_width / window_height), 0.1, 1000.0)
    glMatrixMode(GL_MODELVIEW)

    glPointSize(2.0)

# сохранение в bmp
def save_bmp(filename="output.bmp"):
    try:
        viewport = glGetIntegerv(GL_VIEWPORT)
        x, y, width, height = viewport

        # настройка параметров чтения
        glReadBuffer(GL_BACK)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)

        # чтение пикселей
        pixel = glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE)
        if pixel is None:
            print("Не удалось прочитать данные")
            return False

        # создание изображения
        image = Image.frombytes("RGB", (width, height), pixel)
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        # сохранение
        image.save(filename, "BMP")
        print(f"Изображение сохранено как {filename}")
        return True

    except Exception as e:
        print(f"Ошибка при сохранении BMP: {e}")
        return False


def render(depthmap_array, light, camera, shading, filename="render.bmp"):
    print("Cоздание bmp")
    try:
        # инициализация GLUT
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)

        # создание окна
        glutInitWindowSize(window_width, window_height)
        glutInitWindowPosition(0, 0)
        window_id = glutCreateWindow(b"Hidden Render")

        # настройка OpenGL
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.0, 0.0, 0.0, 0.0)

        glMatrixMode(GL_PROJECTION)
        gluPerspective(45, (window_width / window_height), 0.1, 1000.0)
        glMatrixMode(GL_MODELVIEW)

        # рендер сцены
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        gluLookAt(camera[0], camera[1], camera[2],
                  0, 0, 0,
                  0, 1, 0)
        height, width = depthmap_array.shape
        max_depth = np.max(depthmap_array) if np.max(depthmap_array) != 0 else 1
        scale_factor = 200

        glBegin(GL_QUADS)
        for i in range(height - 1):
            for j in range(width - 1):
                if (depthmap_array[i, j] != 0 and
                        depthmap_array[i, j + 1] != 0 and
                        depthmap_array[i + 1, j] != 0 and
                        depthmap_array[i + 1, j + 1] != 0):

                    vertices = [
                        ((j - width / 2) / width, (height / 2 - i) / height, -depthmap_array[i, j] / max_depth),
                        (((j + 1) - width / 2) / width, (height / 2 - i) / height,
                         -depthmap_array[i, j + 1] / max_depth),
                        (((j + 1) - width / 2) / width, (height / 2 - (i + 1)) / height,
                         -depthmap_array[i + 1, j + 1] / max_depth),
                        ((j - width / 2) / width, (height / 2 - (i + 1)) / height,
                         -depthmap_array[i + 1, j] / max_depth)
                    ]

                    normal = calculate_normal(vertices[0], vertices[1], vertices[2])
                    color = models(normal, light, camera, shading)

                    glColor3f(*color)
                    for vertex in vertices:
                        glVertex3f(vertex[0] * scale_factor, vertex[1] * scale_factor, vertex[2] * scale_factor)
        glEnd()
        # сохранение в BMP
        success = save_bmp(filename)
        glutDestroyWindow(window_id)
        if success:
            print(f"Изображение создано: {filename}")
        return success

    except Exception as e:
        print(f"Ошибка при создании BMP: {e}")
        return False

def main():
    # чтение json файла
    json_file = 'config.json'
    depthmap_file, config = read_json(json_file)
    depthmap_array = read_DeapthMap(depthmap_file)
    light = config['light_source']['position']
    camera = config['viewer']['position']
    model = config.get('reflection_model', 'phong')

    # экспорт
    export_format = config.get('export_format', 'vrml')
    output = config.get('export_name', 'output')
    export(depthmap_array, export_format, output)

    # bmp изображение
    bmp_filename = f"{model}.bmp"
    render(depthmap_array, light, camera, model, bmp_filename)

    # визуализация
    init_glut(
        lambda: render_scene(
            depthmap_array,
            light,
            camera,
            model
        )
    )

    glutMainLoop()

if __name__ == "__main__":
    main()
