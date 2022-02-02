setlocal

@rem enter this directory
cd /d %~dp0

REM integrate with CMake
set PROTOC=..\..\..\vcpkg_5568f11\installed\x64-windows\tools\protobuf\protoc.exe
set PLUGIN=..\..\..\vcpkg_5568f11\installed\x64-windows\tools\grpc\grpc_cpp_plugin.exe

%PROTOC% --cpp_out=. core.proto

endlocal

PAUSE