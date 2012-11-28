!===================
module data_diag
!===================

  implicit none
  public
  logical               :: dump_mesh_info_files
  logical               :: dump_mesh_info_screen
  logical               :: dump_mesh_vtk
  character(len=200)    :: diagpath
  integer               :: lfdiag

!===================
end module data_diag
!===================
