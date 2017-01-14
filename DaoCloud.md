docker-machine satart default
docker-machine ssh default
sudo sed -i "s|EXTRA_ARGS='|EXTRA_ARGS='--registry-mirror=http://7a27ee44.m.daocloud.io |g"/var/lib/boot2docker/profile
exit
docker-machine restart default
