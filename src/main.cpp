#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <boost/algorithm/string.hpp>

#include <igl/pathinfo.h>

#include "Dog/Dog.h"
#include "ModelState.h"
#include "ModelViewer.h"

using namespace std;

bool is_optimizing = true;

ModelState state;
ModelViewer modelViewer(state, state.DC);

const int DEFAULT_GRID_RES = 21;
int editedSubmeshI = -1; // -1 means the entire mesh, i means the i connected component submesh 

void clear_all_and_set_default_params(igl::opengl::glfw::Viewer& viewer) {
  state.DC.init_from_new_dog(state.dog);
}

void run_optimization() {
  state.DC.single_optimization();
}

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch (key) {
  case ' ':
    is_optimizing = !is_optimizing;
    break;
  case 'S':
    state.DC.edit_mode = DogEditor::SELECT_POSITIONAL;
    break;
  case 'D':
    state.DC.edit_mode = DogEditor::TRANSLATE;
    break;
  case 'C':
    state.DC.reset_constraints();
    break;
  case 'R':
    state.DC.single_optimization();
    break;
  case 'E':
    exit(1);
    break;
  }
  return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if ((modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return state.DC.dogEditor->callback_mouse_down();
  return false;
}
bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) {
  if ((modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return  state.DC.dogEditor->callback_mouse_move(mouse_x, mouse_y);
  return false;
}
bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
  if ((modelViewer.viewMode == ViewRulings) || (modelViewer.viewMode == ViewModeMeshWire)) return  state.DC.dogEditor->callback_mouse_up();
  return false;
}

bool callback_pre_draw(igl::opengl::glfw::Viewer& viewer) {
  if ( ((state.DC.has_constraints())) && is_optimizing) run_optimization();
  modelViewer.render(viewer);
  return false;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: dog_editor INPUT_PATH.OBJ or dog_editor planar (optional width and height)" << endl;
    exit(1);
  }
  const std::string input_path = argv[1];
  std::string dirname,basename,extension,filename; igl::pathinfo(input_path, dirname, basename, extension, filename);
  if (boost::iequals(basename, "planar")) {
    int x_res,y_res; x_res = y_res = DEFAULT_GRID_RES;
    if (argc > 2) {x_res = y_res = std::stoi(argv[2]);};
    state.init_from_planar(x_res,y_res);
  } else {
    // Assume obj/off or other types
    state.init_from_mesh(input_path);
  }
  // Set up viewer
  igl::opengl::glfw::Viewer viewer;
  state.DC.init_viewer(viewer);
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
    // Define next window position + size
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(210, 350), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(
        "DOG", nullptr,
        ImGuiWindowFlags_NoSavedSettings
    );
    ImGui::Combo("View mode", (int *)(&modelViewer.viewMode), "MeshWire\0Gauss Map\0Rulings\0\0");
    ImGui::Combo("Edit mode", (int *)(&state.DC.edit_mode), "Select\0Translate\0None\0\0");
    
    ImGui::InputDouble("Bending", &state.DC.p.bending_weight, 0, 0, "%.4f");
    ImGui::InputDouble("Isometry", &state.DC.p.isometry_weight, 0, 0, "%.4f");
    ImGui::InputDouble("Edge regularization", &state.DC.p.reg_edge_weight, 0, 0, "%.4f");
    ImGui::InputDouble("Soft constraints", &state.DC.p.soft_pos_weight, 0, 0, "%.4f");

    ImGui::InputDouble("Constraints deviation", &state.DC.constraints_deviation);
    ImGui::InputDouble("objective", &state.DC.objective);
    ImGui::Checkbox("Is optimizing?", &is_optimizing);

    ImGui::InputDouble("Rulings length", &modelViewer.rulings_length);
    ImGui::InputInt("Rulings modulo", &modelViewer.rulings_mod);

    ImGui::End();
  };
  clear_all_and_set_default_params(viewer);
  
  viewer.data().set_mesh(state.dog.getVrendering(), state.dog.getFrendering());
  viewer.core.align_camera_center(state.dog.getVrendering(), state.dog.getFrendering());
  
  viewer.callback_key_down = callback_key_down;
  viewer.callback_pre_draw = callback_pre_draw; // calls at each frame
  viewer.core.is_animating = true;
  viewer.core.animation_max_fps = 30;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.data().line_width = 2;

  viewer.data().show_lines = false;
  viewer.launch();
}