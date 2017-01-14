docker-machine satart default
docker-machine ssh default
sudo sed -i "s|EXTRA_ARGS='|EXTRA_ARGS='--registry-mirror= |g"/var/lib/boot2docker/profile
exit
docker-machine restart default
