// Delaunay triangulation and Curve construction
// g++ -std=c++14  -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 de_t.cpp -o test

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

using namespace std; 

typedef struct point {
	int x;
	int y;
  unsigned long int id;
}v;

struct nearest_neighbor {
  point p;
  float dist;
};

struct edge{
  point p1, p2;
};

struct pb{
  float slope;
  float c;
};

typedef std::vector<v> vertices; 

struct inner_edge {
  point p1;
  point p2;
  vertices outer_points;
  point outer_point;
  bool visited;
  string id;
  float dist;
};

struct graph_edge{
  point p1, p2;
};

struct triangle{
  point p0, p1, p2;
  int id;
}T;

struct circle{
  float x;
  float y;
  float radius;
}c;

typedef std::vector<inner_edge> edges;
typedef std::vector<edge> c_edges;
typedef std::vector<graph_edge> g_edges;

typedef std::vector<triangle> triangles;
vertices outer_points;

int plot_triangles_and_curve(vertices *points, triangles *tri, c_edges *curve_edges, int argc, char** argv);
void reconstruct_inner_edges(edges *inner_edges, inner_edge e);
void unmark_previous(edges &inner_edges, edges *aa, inner_edge e);
void construct_curve(edges &d_edges, vertices &vert, c_edges &curve_edges);
float distance(v p1, v p2);
float point_orientation(v p1, v p2, v p3){ 
  // >0 for P3 left of the line through P1 and P2
  // =0 for P3 on the line
  // <0 for P3 right of the line
  return (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
}

bool point_inside_circle(v pt, circle c){
  float dist = sqrt((c.x - pt.x) * (c.x - pt.x) +
                        (c.y - pt.y) * (c.y - pt.y));
  return dist < c.radius;
}

string get_id(v p1, v p2){
  char buffer [50];
  int n = sprintf (buffer, "%lu%lu", p1.id, p2.id);
  return buffer;
}

void swap_p(v *p1,v *p2){
  v temp = *p1;
  *p1 = *p2;
  *p2 = temp;
}

void sortx(vertices &points){
  int i, j, n = points.size();
  for (i=0; i< n-1;i++){
    //last i element are already in place
    for(j= 0; j< n-i-1; j++){
      if(points[j].x > points[j+1].x)
        swap_p(&points[j],&points[j+1]);
    }
  } 
}

void sorty(vertices &points){
  int i, j, n = points.size();

  for (i=0; i< n-1;i++){
    //last i element are already in place
    for(j= 0; j< n-i-1; j++){
      if(points[j].y > points[j+1].y)
        swap_p(&points[j],&points[j+1]);
    }
  } 
}

inner_edge create_edge(v p1, v p2, v outer_point, bool visited=false){
  inner_edge e;
  e.visited = visited;
  e.id = get_id(p1, p2);
  e.p1 = p1;  e.p2 = p2; e.outer_point = outer_point;
  return e;
}

circle triangle_circle(v p0, v p1, v p2){
  float ax = p1.x - p0.x;
  float ay = p1.y - p0.y;
  float bx = p2.x - p0.x;
  float by = p2.y - p0.y;

  float m = p1.x * p1.x - p0.x * p0.x + p1.y * p1.y - p0.y * p0.y;
  float u = p2.x * p2.x - p0.x * p0.x + p2.y * p2.y - p0.y * p0.y;
  float s = 1. / (2. * (ax * by - ay * bx));

  circle c;
  c.x = ((p2.y - p0.y) * m + (p0.y - p1.y) * u) * s;
  c.y = ((p0.x - p2.x) * m + (p1.x - p0.x) * u) * s;

  float dx = p0.x - c.x;
  float dy = p0.y - c.y;
  c.radius = sqrt(dx * dx + dy * dy);
  return c;
}

void common_edges(edges &inner_edges, edges &tri_common_edges, triangles &tri){
  unsigned long int id = 1;
  for(int i=0; i < inner_edges.size(); i++){
    for(int j= 0; j < inner_edges.size(); j++ ){
      if( i == j){
        continue;
      }
      if(inner_edges[j].visited == true){
        continue;
      }
      if(inner_edges[i].id == inner_edges[j].id || inner_edges[i].id == get_id(inner_edges[j].p2, inner_edges[j].p1)){
        if(inner_edges[i].outer_point.id == inner_edges[j].outer_point.id){
          continue;
        }
        inner_edges[j].visited = true;
        inner_edge e;
        e.p1 = inner_edges[i].p1;e.p2 = inner_edges[i].p2;
        e.visited = false;
        e.outer_points.push_back(inner_edges[i].outer_point);
        e.outer_points.push_back(inner_edges[j].outer_point);
        tri_common_edges.push_back(e);
        triangle ch_t_a = {e.p1, e.p2, inner_edges[i].outer_point};
        tri.push_back(ch_t_a);
        triangle ch_t_b = {e.p1, e.p2, inner_edges[j].outer_point};
        tri.push_back(ch_t_b);
        break;
      }
    }
    inner_edges[i].visited = true;
  }
}

float distance(v p1, v p2){
  return sqrt((p2.x - p1.x) * (p2.x - p1.x) +
                        (p2.y - p1.y) * (p2.y - p1.y));
}

void delaunay_triangle(edges &inner_edges, triangles &t){
  triangle tr;
  vertices outer_points;
  inner_edge e;
  edges aa;
  while(!inner_edges.empty()){
    e = inner_edges.back();
    inner_edges.pop_back();
    outer_points = e.outer_points;
    point first_point = outer_points[0];
    point second_point = outer_points[1];
    circle c1 = triangle_circle(e.p1, e.p2, first_point);
    circle c2 = triangle_circle(e.p1, e.p2, second_point);

    bool flip = point_inside_circle(second_point, c1) && point_inside_circle(first_point, c2);

    if(flip){
      inner_edge a = create_edge(first_point, second_point, e.p1);
      a.outer_points.push_back(e.p1);
      a.outer_points.push_back(e.p2);
      reconstruct_inner_edges(&inner_edges, a);
      inner_edges.push_back(a);
      unmark_previous(inner_edges, &aa, a);
    }else{
      e.visited = true;
      aa.push_back(e);
    }
  }
  inner_edges.clear();

  for (std::vector<inner_edge>::iterator itr = aa.begin(); itr != aa.end(); ++itr){
    inner_edge in;
    in.p1 = (*itr).p1; in.p2 = (*itr).p2, in.dist = distance((*itr).p1, (*itr).p2);
    inner_edges.push_back(in);
    in.p1 = (*itr).p2; in.p2 = (*itr).outer_points[0], in.dist = distance((*itr).p2, (*itr).outer_points[0]);
    inner_edges.push_back(in);
    in.p1 = (*itr).outer_points[0]; in.p2 = (*itr).p1, in.dist = distance((*itr).outer_points[0], (*itr).p1);
    inner_edges.push_back(in);

    in.p1 = (*itr).p2; in.p2 = (*itr).outer_points[1], in.dist = distance((*itr).p2, (*itr).outer_points[1]);
    inner_edges.push_back(in);
    in.p1 = (*itr).outer_points[1]; in.p2 = (*itr).p1, in.dist = distance((*itr).outer_points[1], (*itr).p1);
    inner_edges.push_back(in);

    tr.p0 = (*itr).p1; tr.p1 = (*itr).p2; 
    tr.p2 = (*itr).outer_points[0]; t.push_back(tr);
    tr.p2 = (*itr).outer_points[1];; t.push_back(tr);
  }
}

void unmark_previous(edges &inner_edges, edges *aa, inner_edge e){
  edges tri_common_edges;
  v flip_o_p1 = e.outer_points[0];
  v flip_o_p2 = e.outer_points[1];
  v flip_p1 = e.p1;
  v flip_p2 = e.p2;
  bool replace = false;
  v replace_p1, replace_p2;
  std::vector<int> index;
  int i = 0;
  for(std::vector<inner_edge>::iterator it = (*aa).begin(); it != (*aa).end(); ++it){
    replace = false;
    v o_p1 = (*it).outer_points[0];
    v o_p2 = (*it).outer_points[1];
    v p1 = (*it).p1;
    v p2 = (*it).p2;

    if(flip_p1.id == p1.id || flip_p1.id == p2.id){
      replace_p1 = flip_p2;
      replace = true;
    }else if(flip_p2.id == p1.id || flip_p2.id == p2.id){
      replace_p1 = flip_p1;
      replace = true;
    }
    if(replace){
      replace = false;
      if(flip_o_p1.id == p1.id || flip_o_p1.id == p2.id){
        replace_p2 = flip_o_p2;
        replace = true;
      }else if(flip_o_p2.id == p1.id || flip_o_p2.id == p2.id){
        replace_p2 = flip_o_p1;
        replace = true;
      }
    }

    if(replace){
      (*it).outer_points.clear();
      (*it).outer_points.push_back(replace_p1);
      if(o_p1.id != replace_p2.id){
        (*it).outer_points.push_back(o_p1);
      }else if(o_p2.id != replace_p2.id){
        (*it).outer_points.push_back(o_p2);
      }
      e.visited = false;
      inner_edges.push_back(*it);
    }else{
      index.push_back(i);
    }
    i++;
  }

  edges temp;
  for(int j=0; j < index.size(); j++){
    temp.push_back((*aa)[index[j]]);
  }
  (*aa).clear();
  *aa = temp;
}

void reconstruct_inner_edges(edges *inner_edges, inner_edge e){
  edges tri_common_edges;
  v flip_o_p1 = e.outer_points[0];
  v flip_o_p2 = e.outer_points[1];
  v flip_p1 = e.p1;
  v flip_p2 = e.p2;
  bool replace = false;
  v replace_p1, replace_p2;

  for(std::vector<inner_edge>::iterator it = (*inner_edges).begin(); it != (*inner_edges).end(); ++it){
    v o_p1 = (*it).outer_points[0];
    v o_p2 = (*it).outer_points[1];
    v p1 = (*it).p1;
    v p2 = (*it).p2;
    replace = false;

    if(flip_p1.id == p1.id || flip_p1.id == p2.id){
      replace_p1 = flip_p2;
      replace = true;
    }else if(flip_p2.id == p1.id || flip_p2.id == p2.id){
      replace_p1 = flip_p1;
      replace = true;
    }
    if(replace){
      replace = false;
      if(flip_o_p1.id == p1.id || flip_o_p1.id == p2.id){
        replace_p2 = flip_o_p2;
        replace = true;
      }else if(flip_o_p2.id == p1.id || flip_o_p2.id == p2.id){
        replace_p2 = flip_o_p1;
        replace = true;
      }
    }

    if(replace){
      (*it).outer_points.clear();
      (*it).outer_points.push_back(replace_p1);
      if(o_p1.id != replace_p2.id){
        (*it).outer_points.push_back(o_p1);
      }else if(o_p2.id != replace_p2.id){
        (*it).outer_points.push_back(o_p2);
      }
      (*it).visited = false;
    }
  }
}

void triangulation(vertices *vertex, triangles *tri, c_edges *curve){
  edges inner_edges;
  unsigned long n;
  sorty(*vertex);
  sortx(*vertex);
  
  vertices upper_chain;
  vertices lower_chain;
  triangles T;
  triangle ch_t;
  float orientation = point_orientation(vertex->at(0), vertex->at(1), vertex->at(2));

  if(orientation != 0){
    ch_t.p0=(*vertex)[0];ch_t.p1=(*vertex)[1];ch_t.p2=(*vertex)[2];
    (*tri).push_back(ch_t);
    inner_edges.push_back(create_edge((*vertex)[0], (*vertex)[1], (*vertex)[2]));
    inner_edges.push_back(create_edge((*vertex)[2], (*vertex)[0], (*vertex)[1]));
    inner_edges.push_back(create_edge((*vertex)[1], (*vertex)[2], (*vertex)[0]));
  }
  
  if(orientation > 0) {
    lower_chain.push_back(vertex->at(0));
    lower_chain.push_back(vertex->at(1));
    lower_chain.push_back(vertex->at(2));

    upper_chain.push_back(vertex->at(0));
    upper_chain.push_back(vertex->at(2));
  }else if(orientation < 0){
    upper_chain.push_back(vertex->at(0));
    upper_chain.push_back(vertex->at(1));
    upper_chain.push_back(vertex->at(2));

    lower_chain.push_back(vertex->at(0));
    lower_chain.push_back(vertex->at(2));
  }else{
    //Collinear
    upper_chain.push_back(vertex->at(0));
    upper_chain.push_back(vertex->at(1));
    upper_chain.push_back(vertex->at(2));

    lower_chain.push_back(vertex->at(0));
    lower_chain.push_back(vertex->at(1));
    lower_chain.push_back(vertex->at(2));
  }

  for (int i=3; i < (*vertex).size(); i++) { 
    v new_point = vertex->at(i);
    v last_point = vertex->at(i-1);

    // Repeated points 
    if(last_point.x == vertex->at(i).x && last_point.y == vertex->at(i).y)
      continue;

    n = upper_chain.size() - 1;
    while(n > 0 && point_orientation(new_point, upper_chain[n], upper_chain[n-1]) < 0) {
      inner_edges.push_back(create_edge(new_point, upper_chain[n], upper_chain[n-1]));
      inner_edges.push_back(create_edge(upper_chain[n-1], new_point, upper_chain[n]));
      inner_edges.push_back(create_edge(upper_chain[n], upper_chain[n-1], new_point));
      ch_t.p0 = new_point; ch_t.p1 = upper_chain[n]; ch_t.p2 = upper_chain[n-1];
      (*tri).push_back(ch_t);
      upper_chain.pop_back();
      n--;
    }

    n = lower_chain.size() - 1;
    while(n > 0 && point_orientation(new_point, lower_chain[n], lower_chain[n-1]) > 0) {
      inner_edges.push_back(create_edge(new_point, lower_chain[n], lower_chain[n-1]));
      inner_edges.push_back(create_edge(lower_chain[n-1], new_point, lower_chain[n]));
      inner_edges.push_back(create_edge(lower_chain[n], lower_chain[n-1], new_point));
      ch_t.p0 = new_point; ch_t.p1 = lower_chain[n]; ch_t.p2 = lower_chain[n-1];
      (*tri).push_back(ch_t);
      lower_chain.pop_back();
      n--;
    }
    upper_chain.push_back(vertex->at(i));
    lower_chain.push_back(vertex->at(i));
  }
  (*tri).clear();
  edges tri_common_edges;
  common_edges(inner_edges, tri_common_edges, *tri);
  (*tri).clear();
  delaunay_triangle(tri_common_edges, *tri);
  construct_curve(tri_common_edges, *vertex, *curve);
}

pb perpendicular_bisector(v p1, v p2){
  pb p;
  float x = 1./2. * (p1.x + p2.x); 
  float y = 1./2. * (p1.y + p2.y);
  float d = (p2.y - p1.y);
  float slope = d /(p2.x - p1.x);
  p.slope = -1. * (1./slope);
  p.c = y - (p.slope * x);
  return p;
}

int orientation_p(pb pbb, v p){
  float val = pbb.slope * p.x + pbb.c;
  if(pbb.slope > 0) {
    if(p.y > val)
      return 1;
    else if(p.y < val)
      return -1;
    else
      return 0;
  } else if (pbb.slope < 0){
    if(p.y < val)
      return 1;
    if(p.y > val)
      return -1;
    else
      return 0;
  }else{
    return 0;
  }
}

void construct_curve(edges &edg, vertices &vert, c_edges &curve_edges){
  for(std::vector<v>::iterator itr = vert.begin(); itr != vert.end(); ++itr){
    nearest_neighbor n;
    std::vector<nearest_neighbor> n_points;
    bool first_point = true;
    n_points.clear();
    for(std::vector<inner_edge>::iterator t_itr = edg.begin(); t_itr != edg.end(); ++t_itr){
      if((*itr).id == (*t_itr).p1.id || (*itr).id == (*t_itr).p2.id){
        v np = (*itr).id == (*t_itr).p1.id ? (*t_itr).p2 : (*t_itr).p1;
        if(first_point){
          n.p = np;
          n.dist = (*t_itr).dist;
          first_point = false;
        }else if((*t_itr).dist < n.dist){ 
          n_points.push_back(n);
          n.p = np;
          n.dist = (*t_itr).dist;
        }else{
          nearest_neighbor npp = {np, (*t_itr).dist};
          n_points.push_back(npp);
        }
      }
    }
    edge e = {*itr, n.p};
    curve_edges.push_back(e);

    bool pb_present = false;
    first_point = true;
    pb p;
    int p_orientation, x;
    float midpoint;
    if(((*itr).x - n.p.x) == 0){
      midpoint = ((*itr).y + n.p.y)/2;
      x = 1;
    }else if(((*itr).y - n.p.y) == 0){
      midpoint = ((*itr).x + n.p.x)/2;
      x = 0;
    }else{
      p = perpendicular_bisector((*itr), n.p);
      pb_present = true;
      p_orientation = orientation_p(p, (*itr));
    }

    for(std::vector<nearest_neighbor>::iterator it = n_points.begin(); it != n_points.end(); ++it){
      bool p_aligned = false;
      if(!pb_present){
        if(x){
          if(((*it).p.y <= midpoint) == ((*itr).y <= midpoint)){ 
            p_aligned = true;
          }
        }else{
          if(((*it).p.x <= midpoint) == ((*itr).x <= midpoint)){ 
            p_aligned = true;
          }
        }
      }else if(orientation_p(p, (*it).p) == p_orientation){
        p_aligned = true;
      }
      if(p_aligned && (*it).p.id != n.p.id){
        if(first_point){
          n = *it;
          first_point = false;
        }else if((*it).dist < n.dist){ 
          n = *it;
        }
      }
    }
    e.p2 = n.p;
    curve_edges.push_back(e);
  }
}

int main(int argc, char** argv) 
{ 
  if(argc != 2){
    printf("Input file not passed");
    return 1;
  }
  FILE * pFile;
  int x, y;
  const char *filename = argv[1];
  pFile = fopen (filename, "r");
  vertices points;

  if (!pFile)
  {
    std::perror("Failed");
    EXIT_FAILURE;
  }
  unsigned long int i = 1;
  while (fscanf(pFile, "P (%d, %d)\n", &x, &y) == 2)
  {
  	v p = {x, y, i++};
    points.push_back(p);
  }

  fclose(pFile);
  triangles tri;
  cout << points.size() << endl;

  if(points.empty()){
    cout << "Empty points" << endl;
    return 1;
  }

  c_edges curve_edges;
  triangulation(&points, &tri, &curve_edges);
  plot_triangles_and_curve(&points, &tri, &curve_edges, argc, argv);        
  return 0;
}

int plot_triangles_and_curve(vertices *points, triangles *tri, c_edges *curve_edges, int argc, char** argv){
  Display *display;
  char *display_name = NULL;
  if((display = XOpenDisplay(display_name)) == NULL) {
      cout << "Could not open display.\n";
      return -1;
  }
  cout << "Connected to X Server " << XDisplayName(display_name) << endl;

  unsigned int screen_number = DefaultScreen(display);
  unsigned int display_width, display_height;
  display_width = DisplayWidth(display, screen_number);
  display_height = DisplayHeight(display, screen_number);

  Window window;
  unsigned int window_width, window_height;
  window_width = (int) display_width / 1.1;
  window_height = (int) display_height / 1.1;
  window = XCreateSimpleWindow(display, RootWindow(display, screen_number), 0, 0, window_width, window_height, 10,
                              BlackPixel(display, screen_number), WhitePixel(display, screen_number));

  Colormap color_map = XDefaultColormap(display, screen_number);
  cout << "Width: " << display_width << ", height: " << display_height << ", screen number: " << screen_number << endl;

  XSizeHints *size_hints = XAllocSizeHints();
  XWMHints *wm_hints = XAllocWMHints();
  XClassHint *class_hint = XAllocClassHint();
  if(size_hints == NULL || wm_hints == NULL || class_hint == NULL) {
      cout << "Some failure with the hints" << endl;
      return 0;
  }

  size_hints -> flags = PPosition | PSize | PMinSize;
  size_hints -> min_height = 60;
  size_hints -> min_width = 60;

  char *win_name_string = (char *) "Example Window";
  char *icon_name_string = (char *) "Icon for example window";

  XTextProperty win_name;
  XTextProperty icon_name;

  XStringListToTextProperty(&win_name_string, 1, &win_name);
  XStringListToTextProperty(&icon_name_string, 1, &icon_name);

  wm_hints -> flags = StateHint | InputHint;
  wm_hints -> initial_state = NormalState;
  wm_hints -> input = False;

  class_hint -> res_name = (char *) "x_use_example";
  class_hint -> res_class = (char *) "examples";

  XSetWMProperties(display, window, &win_name, &icon_name, argv, argc, size_hints, wm_hints, class_hint);
  XSelectInput(display, window, ExposureMask | StructureNotifyMask | ButtonPressMask);

  // Display a window on the screen.
  XMapWindow(display, window);

  XColor temporary_color_1, temporary_color_2;
  XGCValues gc_values;
  GC gc_red = XCreateGC(display, window, 0, &gc_values);
  XSetLineAttributes(display, gc_red, 1, LineSolid, CapRound, JoinRound);
  if(XAllocNamedColor(display, color_map, "red", &temporary_color_1, &temporary_color_2) == 0) {
      cout << "Failed to get color\n";
      return 0;
  }
  else {
      XSetForeground(display, gc_red, temporary_color_1.pixel);
  }

  XColor temporary_color_3, temporary_color_4;
  GC gc_black = XCreateGC(display, window, 0, &gc_values);
  XSetLineAttributes(display, gc_black, 1, LineSolid, CapRound, JoinRound);
  if(XAllocNamedColor(display, color_map, "black", &temporary_color_3, &temporary_color_4) == 0) {
      cout << "Failed to get color\n";
      return 0;
  }
  else {
      XSetForeground(display, gc_black, temporary_color_3.pixel);
  }

  XColor temporary_color_5, temporary_color_6;
  GC gc_blue = XCreateGC(display, window, 0, &gc_values);
  XSetLineAttributes(display, gc_blue, 1, LineSolid, CapRound, JoinRound);
  if(XAllocNamedColor(display, color_map, "blue", &temporary_color_5, &temporary_color_6) == 0) {
      cout << "Failed to get color\n";
      return 0;
  }
  else {
      XSetForeground(display, gc_blue, temporary_color_5.pixel);
  }

  XColor temporary_color_7, temporary_color_8;
  GC gc_yellow = XCreateGC(display, window, 0, &gc_values);
  XSetLineAttributes(display, gc_yellow, 1, LineSolid, CapRound, JoinRound);
  if(XAllocNamedColor(display, color_map, "yellow", &temporary_color_7, &temporary_color_8) == 0) {
      cout << "Failed to get color\n";
      return 0;
  }
  else {
      XSetForeground(display, gc_yellow, temporary_color_7.pixel);
  }
  
  XEvent event;
  while(true) {
      XNextEvent(display, &event);
      switch(event.type) {
          case Expose:
              for (std::vector<v>::iterator it = (*points).begin(); it != (*points).end(); ++it){
                XFillArc(display, window, gc_black, (*it).x - 4,
                        (*it).y - 4, 8, 8, 0, 360 * 64);
              }
              for (std::vector<triangle>::iterator it = (*tri).begin(); it != (*tri).end(); ++it){
                XDrawLine(display, window, gc_blue, (*it).p0.x, (*it).p0.y, (*it).p1.x, (*it).p1.y);
                XDrawLine(display, window, gc_blue, (*it).p1.x, (*it).p1.y, (*it).p2.x, (*it).p2.y);
                XDrawLine(display, window, gc_blue, (*it).p2.x, (*it).p2.y, (*it).p0.x, (*it).p0.y);
              }
              for (std::vector<edge>::iterator it = (*curve_edges).begin(); it != (*curve_edges).end(); ++it){
                XDrawLine(display, window, gc_red, (*it).p1.x, (*it).p1.y, (*it).p2.x, (*it).p2.y);
              }
      }
  }
  return 0;
}
