include_directories (${FINDCLUSTER_SOURCE_DIR}/include)
link_directories (${FINDCLUSTER_BINARY_DIR})

add_executable (findcluster.x findcluster.cpp findcluster_cluster.cpp findcluster_path.cpp findcluster_radius.cpp findcluster_box.cpp findcluster_writemeasure.cpp findcluster_helper.cpp findcluster_write.cpp findcluster_obs.cpp)
