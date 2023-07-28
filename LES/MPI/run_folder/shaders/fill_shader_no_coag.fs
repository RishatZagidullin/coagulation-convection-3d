#version 330 core

out vec4 FragColor;

in vec3 ourColor;

vec3 output;

float val = 1e-2;
float color;

void main()
{
    if (ourColor[0] > val)
    {
        color = ourColor[0]/0.03;
        if (color > 1.) color = 1.;
        output = vec3(color, 0.0, 1. - color);
    }    
    else
    {
        output = vec3(1.,1.,1.);
    }
    FragColor = vec4(output, 1.0);
    
}
