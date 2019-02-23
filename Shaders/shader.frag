#version 330 core
struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shininess;
};

struct DirLight {
    int on;
    vec3 direction;
    
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

struct PointLight {
    int on;
    vec3 position;
    
    float constant;
    float linear;
    float quadratic;
    
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

struct SpotLight {
    int on;
    vec3 position;
    vec3 direction;
    //    float cutOff;
    //    float outerCutOff;
    float cutoff;
    float spot_exponent;
    
    float constant;
    float linear;
    float quadratic;
    
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

#define NR_POINT_LIGHTS 1

in vec3 FragPos;
in vec3 Normal;

out vec4 color;

uniform vec3 viewPos;
uniform DirLight dirLight;
uniform DirLight dirLight2;
uniform PointLight pointLights[NR_POINT_LIGHTS];
//uniform PointLight pointLight;

uniform SpotLight spotLight;
uniform Material material;

uniform bool object_in_collision;
uniform bool drawingBB;


// Function prototypes
vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir);
vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir);
vec3 CalcSpotLight(SpotLight light, vec3 normal, vec3 fragPos, vec3 viewDir);

void main()
{
    // Properties
    vec3 norm = normalize(Normal);
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 result = vec3(0.f,0.f,0.f);
    
    
    //    Dirlight
    if (dirLight.on != 0)
        result += CalcDirLight(dirLight, norm, viewDir);
    
    //PointLight(s)
    for(int i = 0; i < NR_POINT_LIGHTS; i++)
        if (pointLights[i].on != 0)
            result += CalcPointLight(pointLights[i], norm, FragPos, viewDir);
    
    //SpotLight
    if (spotLight.on != 0)
        result += CalcSpotLight(spotLight, norm, FragPos, viewDir);
    
    
    color = vec4(result, 1.0);
//    color = vec4(1,1,1,1);
//    color = vec4(dirLight.specular * material.specular * material.shininess, 1.0f);
}

// Calculates the color when using a directional light.
vec3 CalcDirLight(DirLight light, vec3 normal, vec3 viewDir)
{
    vec3 lightDir = normalize(-light.direction);
    // Diffuse shading
    float diff = max(dot(normal, lightDir), 0.0);
    // Specular shading
    vec3 reflectDir = reflect(-lightDir, normal);
    diff *= dot(viewDir, lightDir);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    
    vec3 ambient = light.ambient * material.ambient;
    vec3 diffuse = light.diffuse * diff * material.diffuse;
    vec3 specular = light.specular * spec * material.specular;
    return (ambient + diffuse + specular);
    //    return material.diffuse;
}

// Calculates the color when using a point light.
vec3 CalcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir)
{
    vec3 lightDir = normalize(light.position - fragPos);
    // Diffuse shading
    float diff = max(dot(normal, lightDir), 0.0);
    // Specular shading
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    // Attenuation
    float distance = length(light.position - fragPos);
    float attenuation = 1.0f / (light.constant + light.linear * distance + light.quadratic * (distance * distance));
    // Combine results
    vec3 ambient = light.ambient * material.ambient;
    vec3 diffuse = light.diffuse * diff * material.diffuse;
    vec3 specular = light.specular * spec * material.specular;
    //    ambient *= attenuation; // Assignment writeup says no
    diffuse *= attenuation;
    specular *= attenuation;
    return (ambient + diffuse + specular);
    //    return (normal); // Debugging
}

// Calculates the color when using a spot light.
vec3 CalcSpotLight(SpotLight light, vec3 normal, vec3 fragPos, vec3 viewDir)
{
    vec3 lightDir = normalize(light.position - fragPos);
    // Diffuse shading
    float diff = max(dot(normal, lightDir), 0.0);
    // Specular shading
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    // Attenuation
    float distance = length(light.position - fragPos);
    float attenuation = 1.0f / (light.constant + light.linear * distance + light.quadratic * (distance * distance));
    // Spotlight intensity
    
    float spotlightdot = dot(light.direction, -lightDir);
    float theta = acos(spotlightdot);
    float intensity = 0;
    if (theta <= light.cutoff) intensity = pow(spotlightdot,light.spot_exponent);
    //    if (true) intensity = pow(spotlightdot,light.spot_exponent);
    
    // Combine results
    vec3 ambient = light.ambient * material.ambient;
    vec3 diffuse = light.diffuse * diff * material.diffuse;
    vec3 specular = light.specular * spec * material.specular;
    //    ambient *= attenuation * intensity; // Assignment writeup says no (actually the discussion formula slides)
    diffuse *= attenuation * intensity;
    specular *= attenuation * intensity;
    return (ambient + diffuse + specular);
}
