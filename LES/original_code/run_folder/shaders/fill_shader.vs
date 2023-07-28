#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

uniform mat4 mvp;
uniform vec4 clip;

out vec3 ourColor;

void main()
{
    gl_Position = mvp*vec4(aPos, 1.0);

    gl_ClipDistance[0] = dot(clip, vec4(aPos, 1.0));
    ourColor = aColor;
}
